import numpy as np
import pandas as pd 
import scipy.interpolate

import subprocess
import os

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt


# запускает из питона стороннюю программу 
def run_comm(command, path):
	subprocess.run(command, 
					encoding="utf-8", shell=True, cwd=f"{path}")


def Sigma_g(R, Sigma0, gama, R_in, R_out):

	Sigma = Sigma0 * (R/R_out)**(-gama) \
			* np.exp(-(R/R_out)**(2-gama)) \
			* np.exp(-(R/R_in)**(gama-2)) 
	
	return Sigma


def Tmid_from_R(meta_name, N):
	# file_path = f"Run33_fast/meta_output/{meta_name}/output_data/r_out-{N:04d}.dat"
	file_path = f"../meta_output/{meta_name}/output_data/r_out-{N:04d}.dat"

	with open(file_path) as input_file:
		time = [next(input_file) for _ in range(2)][1].split()[1]
		input_file.readline()
		input_file.readline()
		head = input_file.readline().split()[1:]
	
	df = pd.read_table(file_path, sep=r"\s+", comment='#', names=head)

	df = df[df.j == 50]
	df.drop(df.tail(1).index, inplace=True) # drop last 1 rows

	T = np.array(df["T"])
	R = np.array(df["R_c"]/AU)

	return T, R


def tau(name, R, Tbase, Tmax, f_d, Sigma0, gama, R_in, R_out):

	# df = pd.read_table(f"../input_opacities/opacity_PR_{name}.txt", sep="\s+", comment="#")
	df = pd.read_table(f"../input_opacities/opacity_PR_{name}.txt", sep="\s+", skiprows=6, names=["temperature[K]", "kappa_P[cm2/g]", "kappa_R[cm2/g]", "kappa_S[cm2/g]"])

	y_interp_P = scipy.interpolate.interp1d(df["temperature[K]"], df["kappa_P[cm2/g]"]/f_d)
	y_interp_R = scipy.interpolate.interp1d(df["temperature[K]"], df["kappa_R[cm2/g]"]/f_d)

	kappaP_arr = np.empty(len(Tbase))
	kappaR_arr = np.empty(len(Tbase))

	for i, t in enumerate(Tbase):
		kappaP_arr[i] = float(y_interp_P(t))
		kappaR_arr[i] = float(y_interp_R(t))


	Sigma_gas = Sigma_g(R, Sigma0, gama, R_in, R_out) / 2
	
	tauP = f_d * Sigma_gas * kappaP_arr
	tauR = f_d * Sigma_gas * kappaR_arr


	x = c_v * (1 + 2*tauP + 3/2*tauP*tauR) / kappaP_arr

	time_cool = x / 4 / SBc / Tmax**3


	Lacc = 450 * Lsun
	F0 = Lacc / 4 / np.pi / (R*AU)**2

	time_heat = x * Tbase / 2 / mu / F0


	return time_cool/secINyear, time_heat/secINyear



def set_ax(ax):
	ax.set_xlabel(r"R [au]")
	ax.set_ylabel(r"time [yr]")
	ax.set_xlim(1e-1,2e2)
	ax.set_ylim(1e-5,1e1)
	ax.tick_params(direction='in', which='both', pad=10, width=2, length=5)
	ax.tick_params(which='major', length=10)
	for side in ax.spines.values():
		side.set_linewidth(2)
	ax.loglog()
	# ax.legend(loc="lower center")
	# ax.legend(loc="center")
	ax.grid(linestyle=':')

	# texts = ax.legend(loc="upper left").get_texts() 
	texts = ax.legend(loc="lower center").get_texts() 

	return texts


def add_legend2(ax):
	ax2 = ax.twinx()
	# ax2.axis("off")
	ax2.set_ylabel(r"time [s]")
	ax2.tick_params(direction='in', which='both', pad=10, width=2, length=5)
	ax2.tick_params(which='major', length=10)
	for side in ax2.spines.values():
		side.set_linewidth(2)
	ax2.set_ylim(1e-5*secINyear,1e1*secINyear)
	ax2.semilogy()


	ax2.plot(-1, -1, "-",  color="black", label=r"t$_{\rm heat}$")
	ax2.plot(-1, -1, "--", color="black", label=r"t$_{\rm cool}$")

	ax2.legend(loc="upper right")


c_light = 2.99792458e10
AU = 1.495978707e13

mu = 0.05
SBc = 5.670367e-5
kboltz = 1.38064852e-16
gamma = 7 / 5 
mmw = 2
amu = 1.660539040e-24

Lsun = 3.828e33

c_v = kboltz / (gamma-1) / mmw / amu

secINyear = 31536000


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
name = "tHC_from_R"
run_comm(f"mkdir {name}", pwd_path)
res = "./" + name


# sigma_gas = Sigma_g(R0, 101.11,  1, 1, 11)
# Mgas = 0
# R0 *= AU
# for i, r in enumerate(R0):
# 	Mgas += 2 * np.pi * r * sigma_gas[i] * (r - R0[i-1])
# Mgas0 = 0.0066 * 1.988416e33
# print(Mgas, Mgas0, Mgas/Mgas0)
Sigma0 = 101.11

plt.rcParams.update({'font.size': 20})


for sca in ["no", "yes"]:
# for sca in ["no"]:

	pdf = PdfPages(f"{res}/{name}_amax_sca{sca}.pdf")

	plt.figure(figsize=(12,10))
	ax = plt.subplot(111)

	for i, amax_micron in enumerate(np.sort(np.append(amax_arr_micron, amax0_micron))):
		meta_name = f"p{p0}f{frac10}amax{amax_micron}sca{sca}"

		# if amax_micron == amax0_micron: k = i

		Tbase, R0 = Tmid_from_R(meta_name, 49)
		Tmax, _ = Tmid_from_R(meta_name, 66)
		t_c, t_h = tau(meta_name, R0, Tbase, Tmax, 0.01, Sigma0, 1, 1, 11)

		if amax_micron == amax0_micron:
			ax.plot(R0, t_h, "-", color=f"C{i}", label=rf"$\mathbf{{a_{{max}} = {amax_micron}\ \mu m}}$")
		else: 
			ax.plot(R0, t_h, "-", color=f"C{i}", label=rf"$a_{{\rm max}} = {amax_micron}\ \mu \rm m$")
		# ax.plot(R0, t_h, "-", color=f"C{i}",  label=r"a$_{\rm max} =$"+f"{amax_micron} $\mu$m")
		ax.plot(R0, t_c, "--", color=f"C{i}")


	texts = set_ax(ax)
	# texts[k].set_weight("bold")

	if sca == "yes":
		plt.title(rf"$p = {p0},\ f_{{\rm Si}} = {frac10},$ with scattering")
	elif sca == "no":
		plt.title(rf"$p = {p0},\ f_{{\rm Si}} = {frac10},$ without scattering")
	# plt.title(f"p = {p0}, fracSi = {frac10}, sca = {sca}")
	add_legend2(ax)

	plt.tight_layout()
	pdf.savefig()
	plt.close()

	pdf.close()

# ----------------------------------------------------------------------------------------------
	pdf = PdfPages(f"{res}/{name}_fracSi_sca{sca}.pdf")

	plt.figure(figsize=(12,10))
	ax = plt.subplot(111)

	for i, frac1 in enumerate(np.sort(np.append(frac1_arr, frac10))):
		meta_name = f"p{p0}f{frac1}amax{amax0_micron}sca{sca}"

		# if frac1 == frac10: k = i

		Tbase, R0 = Tmid_from_R(meta_name, 49)
		Tmax, _ = Tmid_from_R(meta_name, 66)
		t_c, t_h = tau(meta_name, R0, Tbase, Tmax, 0.01, Sigma0, 1, 1, 11)

		if frac1 == frac10:
			ax.plot(R0, t_h, "-", color=f"C{i}", label=rf"$\mathbf{{f_{{Si}} = {frac1}}}$")
		else:
			ax.plot(R0, t_h, "-", color=f"C{i}", label=rf"$f_{{\rm Si}} = {frac1}$")
		# ax.plot(R0, t_h, "-", color=f"C{i}", label=r"frac$_{\rm Si} =$"+f"{frac1}")
		ax.plot(R0, t_c, "--", color=f"C{i}")


	texts = set_ax(ax)
	# texts[k].set_weight("bold")

	# plt.title(f"p = {p0}, amax = {amax0_micron} $\mu$m, sca = {sca}")
	if sca == "yes":
		plt.title(rf"$p = {p0},\ a_{{\rm max}} = {amax0_micron}\ \mu \rm m,$ with scattering")
	elif sca == "no":
		plt.title(rf"$p = {p0},\ a_{{\rm max}} = {amax0_micron}\ \mu \rm m,$ without scattering")
	add_legend2(ax)

	plt.tight_layout()
	pdf.savefig()
	plt.close()

	pdf.close()

# ----------------------------------------------------------------------------------------------
	pdf = PdfPages(f"{res}/{name}_p_sca{sca}.pdf")

	plt.figure(figsize=(12,10))
	ax = plt.subplot(111)

	for i, p in enumerate(np.sort(np.append(p_arr, p0))):
		meta_name = f"p{p}f{frac10}amax{amax0_micron}sca{sca}"

		# if p == p0: k = i

		Tbase, R0 = Tmid_from_R(meta_name, 49)
		Tmax, _ = Tmid_from_R(meta_name, 66)
		t_c, t_h = tau(meta_name, R0, Tbase, Tmax, 0.01, Sigma0, 1, 1, 11)

		if p == p0:
			ax.plot(R0, t_h, "-", color=f"C{i}", label=rf"$\mathbf{{p = {p}}}$")
		else:
			ax.plot(R0, t_h, "-", color=f"C{i}", label=rf"$p = {p}$")
		# ax.plot(R0, t_h, "-", color=f"C{i}", label=f"p = {p}")
		ax.plot(R0, t_c, "--", color=f"C{i}")


	texts = set_ax(ax)
	# texts[k].set_weight("bold")

	# plt.title(f"fracSi = {frac10}, amax = {amax0_micron} $\mu$m, sca = {sca}")
	if sca == "yes":
		plt.title(rf"$f_{{\rm Si}} = {frac10},\ a_{{\rm max}} = {amax0_micron}\ \mu \rm m,$ with scattering")
	elif sca == "no":
		plt.title(rf"$f_{{\rm Si}} = {frac10},\ a_{{\rm max}} = {amax0_micron}\ \mu \rm m,$ without scattering")
	add_legend2(ax)

	plt.tight_layout()
	pdf.savefig()
	plt.close()

	pdf.close()


