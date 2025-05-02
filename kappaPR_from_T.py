

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
	

def get_data(name):
	# df = pd.read_table(f"../input_opacities/opacity_PR_{name}.txt", sep="\s+", comment="#")
	df = pd.read_table(f"../input_opacities/opacity_PR_{name}.txt", sep="\s+", skiprows=6, names=["temperature[K]", "kappa_P[cm2/g]", "kappa_R[cm2/g]", "kappa_S[cm2/g]"])

	return np.array(df["temperature[K]"]), np.array(df["kappa_P[cm2/g]"]), np.array(df["kappa_R[cm2/g]"])


def set_ax(ax):
	ax.set_xlabel(r"T [K]")
	ax.set_ylabel(r"$\varkappa_{\rm P}, \varkappa_{\rm R}$ [cm$^2$ g$_{\rm gas}^{-1}$]")
	ax.set_xlim(1e0,1e4)
	ax.set_ylim(1e-3,1e2)
	ax.tick_params(direction='in', which='both', pad=10, width=2, length=5)
	ax.tick_params(which='major', length=10)
	for side in ax.spines.values():
		side.set_linewidth(2)
	ax.loglog()
	# ax.legend(loc="lower center")
	# ax.legend(loc="center")
	ax.grid(linestyle=':')

	texts = ax.legend(loc="upper left").get_texts() 
	# texts = ax.legend(loc="lower center").get_texts() 

	return texts


def add_text(ax):
	bbox = dict(boxstyle="round", fc="white", ec="gray", alpha=0.6)
	txt = r"$-\ \varkappa_{\rm P}$"+"\n"+r"-- $\varkappa_{\rm R}$"
	ax.text(0.03, 0.97, txt, va="top", ha="left", transform=ax.transAxes, alpha=0, bbox=bbox, rasterized=True)
	ax.text(0.03, 0.97, txt, va="top", ha="left", transform=ax.transAxes)

def add_legend2(ax):
	ax2 = ax.twinx()
	ax2.axis("off")

	ax2.plot(-1, -1, "-",  color="black", label=r"$\varkappa_{\rm P}$")
	ax2.plot(-1, -1, "--", color="black", label=r"$\varkappa_{\rm R}$")

	ax2.legend(loc="upper center")
	


def kappa_law(kappaI_arr, a_arr, b_arr, rho):

	xT = [None] * (len(kappaI_arr)-1)
	# xT = [1] + xT + [1e4]

	for i in range(len(kappaI_arr)-1):
		xT[i] = (kappaI_arr[i+1] / kappaI_arr[i] * rho**(a_arr[i+1]-a_arr[i])) ** (1/(b_arr[i]-b_arr[i+1]))

	xT = [1] + xT + [1e4]

	T_arr = np.empty(0)
	kappa_arr = np.empty(0)

	for i in range(len(xT)-1):
		T = [xT[i]] + [t for t in range(int(xT[i]), int(xT[i+1]))] + [xT[i+1]]
		kappa = kappaI_arr[i] * rho**a_arr[i] * np.array(T)**b_arr[i]

		T_arr = np.concatenate((T_arr, T)) 
		kappa_arr = np.concatenate((kappa_arr, kappa)) 

	return T_arr, kappa_arr


def draw_kappa_law(ax, kappaI_arr, alpha_arr, beta_arr):

	ax2 = ax.twinx()
	ax2.set_xlim(1e0,1e4)
	ax2.set_ylim(5e-4,1e2)
	ax2.axis("off")
	
	# for i, rho in enumerate([1e-5, 1e-9]):
	# 	T_arr, kappa_arr = kappa_law(kappaI_arr, alpha_arr, beta_arr, rho)
	# 	ax.plot(T_arr, kappa_arr, ":", lw=3, color="black", label=r"$\rho =$"+f"10$^{{{int(np.log10(rho))}}} \frac{{g}}{{cm}^3}$")	

	T_arr, kappa_arr = kappa_law(kappaI_arr, alpha_arr, beta_arr, 1e-5)
	ax2.plot(T_arr, kappa_arr, linestyle=(0,(1,1)), lw=3, color="black", label=r"$\rho_{\rm gas} = 10^{-5} \frac{\rm g}{\rm cm^3}$")	
	# df = pd.DataFrame(names=[""])

	T_arr, kappa_arr = kappa_law(kappaI_arr, alpha_arr, beta_arr, 1e-9)
	ax2.plot(T_arr, kappa_arr, linestyle=(0,(1,4)), lw=3, color="black", label=r"$\rho_{\rm gas} = 10^{-9} \frac{\rm g}{\rm cm^3}$")	

	ax2.loglog()
	ax2.legend(loc="lower right")
	# ax2.legend(loc="lower center")
	# ax2.legend(loc=(0.25,0.0165))
	


c_light = 2.99792458e10
AU = 1.495978707e13

kappaI_arr = np.array([2e-4, 2e16, 0.1, 2e81, 1e-8, 1e-36, 1.5e20, 0.348])
alpha_arr = np.array([0, 0, 0, 1, 2/3, 1/3, 1, 0])
beta_arr = np.array([2, -7, 1/2, -24, 3, 10, -5/2, 0])


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
name = "kappaPR_from_T"
run_comm(f"mkdir {name}", pwd_path)
res = "./" + name


plt.rcParams.update({'font.size': 30})


for sca in ["no", "yes"]:
# for sca in ["no"]:

	pdf = PdfPages(f"{res}/{name}_amax_sca{sca}.pdf")

	plt.figure(figsize=(12,10))
	ax = plt.subplot(111)

	draw_kappa_law(ax, kappaI_arr, alpha_arr, beta_arr)

	for i, amax_micron in enumerate(np.sort(np.append(amax_arr_micron, amax0_micron))):
		meta_name = f"p{p0}f{frac10}amax{amax_micron}sca{sca}"

		# if amax_micron == amax0_micron: k = i
		if i in [1,2]: continue

		T, P, R = get_data(meta_name)

		if amax_micron == amax0_micron:
			ax.plot(T, P, "-", color=f"C{i}", label=rf"$\mathbf{{a_{{max}} = {amax_micron}\ \mu m}}$")
		else: 
			ax.plot(T, P, "-", color=f"C{i}", label=rf"$a_{{\rm max}} = {amax_micron}\ \mu \rm m$")
		# ax.plot(T, P, "-", color=f"C{i}",  label=r"a$_{\rm max} =$"+f"{amax_micron} $\mu$m")
		ax.plot(T, R, "--", color=f"C{i}")
		
	texts = set_ax(ax)
	# texts[k].set_weight("bold")

	# plt.title(f"p = {p0}, fracSi = {frac10}, sca = {sca}")
	if sca == "yes":
		plt.title(rf"$p = {p0},\ f_{{\rm Si}} = {frac10},$ with scattering")
	elif sca == "no":
		plt.title(rf"$p = {p0},\ f_{{\rm Si}} = {frac10},$ without scattering")
	# add_text(ax)
	add_legend2(ax)

	plt.tight_layout()
	pdf.savefig()
	plt.close()

	pdf.close()

# ----------------------------------------------------------------------------------------------
	pdf = PdfPages(f"{res}/{name}_fracSi_sca{sca}.pdf")

	plt.figure(figsize=(12,10))
	ax = plt.subplot(111)

	draw_kappa_law(ax, kappaI_arr, alpha_arr, beta_arr)

	for i, frac1 in enumerate(np.sort(np.append(frac1_arr, frac10))):
		meta_name = f"p{p0}f{frac1}amax{amax0_micron}sca{sca}"

		# if frac1 == frac10: k = i

		T, P, R = get_data(meta_name)

		if frac1 == frac10:
			ax.plot(T, P, "-", color=f"C{i}", label=rf"$\mathbf{{f_{{Si}} = {frac1}}}$")
		else:
			ax.plot(T, P, "-", color=f"C{i}", label=rf"$f_{{\rm Si}} = {frac1}$")
		# ax.plot(T, P, "-", color=f"C{i}", label=r"frac$_{\rm Si} =$"+f"{frac1}")
		ax.plot(T, R, "--", color=f"C{i}")

	texts = set_ax(ax)
	# texts[k].set_weight("bold")

	# plt.title(f"p = {p0}, amax = {amax0_micron} $\mu$m, sca = {sca}")
	if sca == "yes":
		plt.title(rf"$p = {p0},\ a_{{\rm max}} = {amax0_micron}\ \mu \rm m,$ with scattering")
	elif sca == "no":
		plt.title(rf"$p = {p0},\ a_{{\rm max}} = {amax0_micron}\ \mu \rm m,$ without scattering")
	# add_text(ax)
	add_legend2(ax)

	plt.tight_layout()
	pdf.savefig()
	plt.close()

	pdf.close()

# ----------------------------------------------------------------------------------------------
	pdf = PdfPages(f"{res}/{name}_p_sca{sca}.pdf")

	plt.figure(figsize=(12,10))
	ax = plt.subplot(111)

	draw_kappa_law(ax, kappaI_arr, alpha_arr, beta_arr)

	for i, p in enumerate(np.sort(np.append(p_arr, p0))):
		meta_name = f"p{p}f{frac10}amax{amax0_micron}sca{sca}"

		# if p == p0: k = i

		T, P, R = get_data(meta_name)

		if p == p0:
			ax.plot(T, P, "-", color=f"C{i}", label=rf"$\mathbf{{p = {p}}}$")
		else:
			ax.plot(T, P, "-", color=f"C{i}", label=rf"$p = {p}$")
		# ax.plot(T, P, "-", color=f"C{i}", label=f"p = {p}")
		ax.plot(T, R, "--", color=f"C{i}")

	texts = set_ax(ax)
	# texts[k].set_weight("bold")

	# plt.title(f"fracSi = {frac10}, amax = {amax0_micron} $\mu$m, sca = {sca}")
	if sca == "yes":
		plt.title(rf"$f_{{\rm Si}} = {frac10},\ a_{{\rm max}} = {amax0_micron}\ \mu \rm m,$ with scattering")
	elif sca == "no":
		plt.title(rf"$f_{{\rm Si}} = {frac10},\ a_{{\rm max}} = {amax0_micron}\ \mu \rm m,$ without scattering")
	# add_text(ax)
	add_legend2(ax)

	plt.tight_layout()
	pdf.savefig()
	plt.close()
	
	pdf.close()




# with open("../meta_names0.txt") as f:
# 	meta_names = f.readlines()

# plt.rcParams.update({'font.size': 30})

# pdf = PdfPages(f"kappaPR_from_T.pdf")

# for name in meta_names:
# 	name = name[:-1]

	
# 	plt.figure(figsize=(12,10))	
# 	ax = plt.subplot(111)

# 	df = pd.read_table(f"../input_opacities/opacity_PR_{name}.txt", sep="\s+", comment="#")
# 	# df = pd.read_table(f"../input_opacities/opacity_PR_{name}.txt", sep="\s+", skiprows=6, names=["temperature[K]", "kappa_P[cm2/g]", "kappa_R[cm2/g]", "kappa_S[cm2/g]"])


# 	ax.plot(df["temperature[K]"], df["kappa_P[cm2/g]"], label="$\kappa_{\rm P}$")
# 	ax.plot(df["temperature[K]"], df["kappa_R[cm2/g]"], label="$\kappa_{\rm R}$")

# 	set_ax(ax)

# 	plt.title(name)
	
# 	plt.tight_layout()
# 	pdf.savefig()
# 	plt.close()


# pdf.close()
