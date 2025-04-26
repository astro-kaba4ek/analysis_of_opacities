

import numpy as np
import pandas as pd 

import subprocess
import os

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt


# запускает из питона стороннюю программу 
def run_comm(command, path):
	subprocess.run(command, 
					encoding="utf-8", shell=True, cwd=f"{path}")
	

def get_data(time_slices, df0, lambda0, meta_name):

	times = np.empty(len(time_slices))
	flux = np.empty(len(time_slices))

	dist = 402 
	
	# file_path = f"Run33_fast/meta_output/{meta_name}/output_data/RADMC3D/0049/SED_HURAKAN.dat"
	file_path = f"../meta_output/{meta_name}/output_data/RADMC3D/0049/SED_HURAKAN.dat"

	df = pd.read_table(file_path, skiprows=3, sep="\s+", names=["lambda", "flux"])
	lambd = np.array(df["lambda"])
	sed = np.array(df["flux"]) * c_light/lambd/1e-4 / dist**2
	
	for i, l in enumerate(lambd):
		if (l - lambda0) == 0: 
			flux0 = sed[i]
			break
		elif (l - lambda0) > 0: 
			flux0 = (lambda0-lambd[i-1]) / (lambd[i]-lambd[i-1]) * (sed[i]-sed[i-1]) + sed[i-1]
			break


	for k, N in enumerate(time_slices):

		# file_path = f"Run33_fast/meta_output/{meta_name}/output_data/RADMC3D/{N}/SED_HURAKAN.dat"
		file_path = f"../meta_output/{meta_name}/output_data/RADMC3D/{N}/SED_HURAKAN.dat"

		df = pd.read_table(file_path, skiprows=3, sep="\s+", names=["lambda", "flux"])
		lambd = np.array(df["lambda"])
		sed = np.array(df["flux"]) * c_light/lambd/1e-4 / dist**2
		
		for i, l in enumerate(lambd):
			if (l - lambda0) == 0: 
				y = sed[i]
				break
			elif (l - lambda0) > 0: 
				y = (lambda0-lambd[i-1]) / (lambd[i]-lambd[i-1]) * (sed[i]-sed[i-1]) + sed[i-1]
				break

		flux[k] = y
		if k == 0:
			times[k] = df0["time[yr]"][int(N)]
		else:
			times[k] = df0["time[yr]"][int(N)-1]

	flux /= flux0

	return times, flux


def set_ax(ax):
	shift = 478
	ax.set_xlabel(r"time [yr]")
	ax.set_ylabel(r"$F_\nu (t)\ /\ F_\nu (t_0)$")
	ax.set_xlim(shift-15, shift+155)
	ax.set_ylim(0.7, 350)
	ax.tick_params(direction='in', which='both', pad=10, width=2, length=5)
	ax.tick_params(which='major', length=10)
	ax.tick_params(axis="x", length=7)
	for side in ax.spines.values():
		side.set_linewidth(2)
	ax.semilogy()
	ax.legend(loc="upper right")
	ax.grid(linestyle=':')

	texts = ax.legend().get_texts() 

	return texts

def add_text(ax, lambda0):
	bbox = dict(boxstyle="round", fc="white", ec="gray", alpha=0.6)
	txt = f"$\lambda =${lambda0:5.1f} $\mu$m"
	ax.text(0.5, 0.03, txt, va="bottom", ha="center", transform=ax.transAxes, alpha=0, bbox=bbox, rasterized=True)
	ax.text(0.5, 0.03, txt, va="bottom", ha="center", transform=ax.transAxes)


c_light = 2.99792458e10
AU = 1.495978707e13


with open("../opacities_for_HURAKAN/param_space.txt") as f:
	lines = f.readlines()

for l in lines:
	l_arr = l.split()[:-1]

	if l_arr[0] == "amax":
		amax0_micron = float(l_arr[1])
		amax_arr_micron = np.sort(np.array(list(map(float, l_arr[2:]))))
		amax0 = amax0_micron * 1e-4
		amax_arr = amax_arr_micron * 1e-4
	elif l_arr[0] == "frac1":
		frac10 = float(l_arr[1])
		frac1_arr =  np.sort(np.array(list(map(float, l_arr[2:]))))
	elif l_arr[0] == "p":
		p0 = float(l_arr[1])
		p_arr =  np.sort(np.array(list(map(float, l_arr[2:]))))



pwd_path = f"{os.getcwd()}"
run_comm("mkdir flux_from_time", pwd_path)
res = "./flux_from_time"

lambda0_arr = [4.6, 100, 850]
# lambda0 = 100 # micron

df0 = pd.read_table("../hurakan/flare_profile_FUOri.inp", sep="\s+")


time_slices = sorted(os.listdir(f"../meta_output/p{p0}f{frac10}amax{amax0_micron}scano/output_data/RADMC3D"))
# time_slices = sorted(os.listdir("Run33_fast/meta_output/p2.5f0.8amax1.0scano/output_data/RADMC3D"))
# time_slices = ["0049", "0066", "0080", "0090", "0115", "0150", "0161", "0181", "0200", "0220", "0240", "0260", "0280"]
# time_slices = ["0049"]

for sca in ["no", "yes"]:
# for sca in ["no"]:

	plt.rcParams.update({'font.size': 20})

	pdf = PdfPages(f"{res}/flux_from_time_amax_sca{sca}.pdf")


	for lambda0 in lambda0_arr:
		plt.figure(figsize=(12,10))
		ax = plt.subplot(111)
			
		for i, amax_micron in enumerate(np.sort(np.append(amax_arr_micron, amax0_micron))):
			meta_name = f"p{p0}f{frac10}amax{amax_micron}sca{sca}"

			# if amax_micron == amax0_micron: k = i

			try:
				times, flux = get_data(time_slices, df0, lambda0, meta_name)			
				if amax_micron == amax0_micron:
					ax.plot(times, flux, label=rf"$\mathbf{{a_{{max}} = {amax_micron}\ \mu m}}$")
				else: 
					ax.plot(times, flux, label=rf"$a_{{\rm max}} = {amax_micron}\ \mu \rm m$")
			except:
				print(f"Cannot create pdf: {meta_name}")

		texts = set_ax(ax)
		# texts[k].set_weight("bold")

		if sca == "yes":
			plt.title(rf"$p = {p0},\ f_{{\rm Si}} = {frac10},$ with scattering")
		elif sca == "no":
			plt.title(rf"$p = {p0},\ f_{{\rm Si}} = {frac10},$ without scattering")
		
		add_text(ax, lambda0)

		plt.tight_layout()
		pdf.savefig()
		plt.close()

	pdf.close()

# ----------------------------------------------------------------------------------------------
	pdf = PdfPages(f"{res}/flux_from_time_fracSi_sca{sca}.pdf")

	for lambda0 in lambda0_arr:
		plt.figure(figsize=(12,10))
		ax = plt.subplot(111)
			
		for i, frac1 in enumerate(np.sort(np.append(frac1_arr, frac10))):
			meta_name = f"p{p0}f{frac1}amax{amax0_micron}sca{sca}"

			# if frac1 == frac10: k = i

			try:
				times, flux = get_data(time_slices, df0, lambda0, meta_name)
				if frac1 == frac10:
					ax.plot(times, flux, label=rf"$\mathbf{{f_{{Si}} = {frac1}}}$")
				else:
					ax.plot(times, flux, label=rf"$f_{{\rm Si}} = {frac1}$")
			except:
				print(f"Cannot create pdf: {meta_name}")

		texts = set_ax(ax)
		# texts[k].set_weight("bold")

		# plt.title(f"p = {p0}, amax = {amax0_micron} $\mu$m, sca = {sca}")
		if sca == "yes":
			plt.title(rf"$p = {p0},\ a_{{\rm max}} = {amax0_micron}\ \mu \rm m,$ with scattering")
		elif sca == "no":
			plt.title(rf"$p = {p0},\ a_{{\rm max}} = {amax0_micron}\ \mu \rm m,$ without scattering")
		
		add_text(ax, lambda0)

		plt.tight_layout()
		pdf.savefig()
		plt.close()

	pdf.close()

# ----------------------------------------------------------------------------------------------
	pdf = PdfPages(f"{res}/flux_from_time_p_sca{sca}.pdf")

	for lambda0 in lambda0_arr:
		plt.figure(figsize=(12,10))
		ax = plt.subplot(111)
			
		for i, p in enumerate(np.sort(np.append(p_arr, p0))):
			meta_name = f"p{p}f{frac10}amax{amax0_micron}sca{sca}"

			# if p == p0: k = i

			try:
				times, flux = get_data(time_slices, df0, lambda0, meta_name)		
				if p == p0:
					ax.plot(times, flux, label=rf"$\mathbf{{p = {p}}}$")
				else:
					ax.plot(times, flux, label=rf"$p = {p}$")
			except:
				print(f"Cannot create pdf: {meta_name}")

		texts = set_ax(ax)
		# texts[k].set_weight("bold")

		# plt.title(f"fracSi = {frac10}, amax = {amax0_micron} $\mu$m, sca = {sca}")
		if sca == "yes":
			plt.title(rf"$f_{{\rm Si}} = {frac10},\ a_{{\rm max}} = {amax0_micron}\ \mu \rm m,$ with scattering")
		elif sca == "no":
			plt.title(rf"$f_{{\rm Si}} = {frac10},\ a_{{\rm max}} = {amax0_micron}\ \mu \rm m,$ without scattering")
		
		add_text(ax, lambda0)

		plt.tight_layout()
		pdf.savefig()
		plt.close()

	pdf.close()