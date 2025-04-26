

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
	

def get_data(meta_name, N):
	dist = 402

	# file_path = f"Run33_fast/meta_output/{meta_name}/output_data/RADMC3D/{N}/SED_HURAKAN.dat"
	file_path = f"../meta_output/{meta_name}/output_data/RADMC3D/{N}/SED_HURAKAN.dat"

	df = pd.read_table(file_path, skiprows=3, sep="\s+", names=["lambda", "flux"])
	lambd = np.array(df["lambda"])
	sed = np.array(df["flux"]) * c_light/lambd/1e-4 / dist**2

	return lambd, sed


def set_ax(ax):
	ax.set_xlabel(r"$\lambda$ [$\mu$m]")
	ax.set_ylabel(r"$\nu F_\nu$ [erg s$^{-1}$ cm$^{-2}$]")
	ax.set_xlim(1e-1,1e4)
	ax.set_ylim(1e-14,1e-7)
	ax.tick_params(direction='in', which='both', pad=10, width=2, length=5)
	ax.tick_params(which='major', length=10)
	for side in ax.spines.values():
		side.set_linewidth(2)
	ax.loglog()
	ax.legend(loc="upper right")
	ax.grid(linestyle=':')

	texts = ax.legend().get_texts() 

	return texts

def add_txt(ax, t_slice):
	bbox = dict(boxstyle="round", fc="white", ec="gray", alpha=0.6)
	if int(t_slice) == 0:
		x = df0["time[yr]"][0]
	else:
		x = df0["time[yr]"][int(t_slice)-1]
	txt = f"Time = {x:6.2f} yr\nTime slice: {t_slice}"
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
run_comm("mkdir SED", pwd_path)
res = "./SED"

df0 = pd.read_table("../hurakan/flare_profile_FUOri.inp", sep="\s+")


time_slices = sorted(os.listdir(f"../meta_output/p{p0}f{frac10}amax{amax0_micron}scano/output_data/RADMC3D"))
# time_slices = sorted(os.listdir("Run33_fast/meta_output/p2.5f0.8amax1.0scano/output_data/RADMC3D"))
# time_slices = ["0049", "0066", "0080", "0090", "0115", "0150", "0161", "0181"]
# time_slices = ["0049"]

for sca in ["no", "yes"]:
# for sca in ["no"]:

	plt.rcParams.update({'font.size': 20})

	pdf = PdfPages(f"{res}/SED_amax_sca{sca}.pdf")


	for t_slice in time_slices:
		plt.figure(figsize=(12,10))
		ax = plt.subplot(111)

		for i, amax_micron in enumerate(np.sort(np.append(amax_arr_micron, amax0_micron))):
			meta_name = f"p{p0}f{frac10}amax{amax_micron}sca{sca}"

			# if amax_micron == amax0_micron: k = i

			try:
				lambd, sed = get_data(meta_name, t_slice)		
				if amax_micron == amax0_micron:
					ax.plot(lambd, sed, label=rf"$\mathbf{{a_{{max}} = {amax_micron}\ \mu m}}$")
				else: 
					ax.plot(lambd, sed, label=rf"$a_{{\rm max}} = {amax_micron}\ \mu \rm m$")	
				# ax.plot(lambd, sed, label=r"a$_{\rm max} =$"+f"{amax_micron} $\mu$m")
			except:
				print(f"Cannot create pdf: {meta_name}")

		texts = set_ax(ax)
		# texts[k].set_weight("bold")

		# plt.title(f"p = {p0}, fracSi = {frac10}, sca = {sca}")
		if sca == "yes":
			plt.title(rf"$p = {p0},\ f_{{\rm Si}} = {frac10},$ with scattering")
		elif sca == "no":
			plt.title(rf"$p = {p0},\ f_{{\rm Si}} = {frac10},$ without scattering")
		add_txt(ax, t_slice)

		plt.tight_layout()
		pdf.savefig()
		plt.close()

	pdf.close()

# ----------------------------------------------------------------------------------------------
	pdf = PdfPages(f"{res}/SED_fracSi_sca{sca}.pdf")

	for t_slice in time_slices:
		plt.figure(figsize=(12,10))
		ax = plt.subplot(111)
			
		for i, frac1 in enumerate(np.sort(np.append(frac1_arr, frac10))):
			meta_name = f"p{p0}f{frac1}amax{amax0_micron}sca{sca}"

			# if frac1 == frac10: k = i

			try:
				lambd, sed = get_data(meta_name, t_slice)		
				if frac1 == frac10:
					ax.plot(lambd, sed, label=rf"$\mathbf{{f_{{Si}} = {frac1}}}$")
				else:
					ax.plot(lambd, sed, label=rf"$f_{{\rm Si}} = {frac1}$")	
				# ax.plot(lambd, sed, label=r"frac$_{\rm Si} =$"+f"{frac1}")
			except:
				print(f"Cannot create pdf: {meta_name}")

		texts = set_ax(ax)
		# texts[k].set_weight("bold")

		# plt.title(f"p = {p0}, amax = {amax0_micron} $\mu$m, sca = {sca}")
		if sca == "yes":
			plt.title(rf"$p = {p0},\ a_{{\rm max}} = {amax0_micron}\ \mu \rm m,$ with scattering")
		elif sca == "no":
			plt.title(rf"$p = {p0},\ a_{{\rm max}} = {amax0_micron}\ \mu \rm m,$ without scattering")
		add_txt(ax, t_slice)

		plt.tight_layout()
		pdf.savefig()
		plt.close()

	pdf.close()

# ----------------------------------------------------------------------------------------------
	pdf = PdfPages(f"{res}/SED_p_sca{sca}.pdf")

	for t_slice in time_slices:
		plt.figure(figsize=(12,10))
		ax = plt.subplot(111)
			
		for i, p in enumerate(np.sort(np.append(p_arr, p0))):
			meta_name = f"p{p}f{frac10}amax{amax0_micron}sca{sca}"

			# if p == p0: k = i

			try:
				lambd, sed = get_data(meta_name, t_slice)	
				if p == p0:
					ax.plot(lambd, sed, label=rf"$\mathbf{{p = {p}}}$")
				else:
					ax.plot(lambd, sed, label=rf"$p = {p}$")		
				# ax.plot(lambd, sed, label=f"p = {p}")
			except:
				print(f"Cannot create pdf: {meta_name}")

		texts = set_ax(ax)
		# texts[k].set_weight("bold")

		# plt.title(f"fracSi = {frac10}, amax = {amax0_micron} $\mu$m, sca = {sca}")
		if sca == "yes":
			plt.title(rf"$f_{{\rm Si}} = {frac10},\ a_{{\rm max}} = {amax0_micron}\ \mu \rm m,$ with scattering")
		elif sca == "no":
			plt.title(rf"$f_{{\rm Si}} = {frac10},\ a_{{\rm max}} = {amax0_micron}\ \mu \rm m,$ without scattering")
		add_txt(ax, t_slice)

		plt.tight_layout()
		pdf.savefig()
		plt.close()

	pdf.close()