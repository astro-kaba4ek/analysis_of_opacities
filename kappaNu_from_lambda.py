
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
	# df = pd.read_table(f"../input_opacities/opacity_nu_{name}.txt", sep="\s+", comment="#")
	df = pd.read_table(f"../input_opacities/opacity_nu_{name}.txt", sep="\s+", skiprows=6, names=["freq[1/s]", "kappa(abs)[cm2/g]", "sigma(sca)[cm2/g]"])

	return c_light/np.array(df["freq[1/s]"]), np.array(df["kappa(abs)[cm2/g]"]), np.array(df["sigma(sca)[cm2/g]"])


def set_ax(ax):
	ax.set_xlabel(r"$\lambda$ [cm]")
	ax.set_ylabel(r"$\varkappa_\nu^{\rm abs}, \varkappa_\nu^{\rm sca}$ [cm$^2$ g$_{\rm dust}^{-1}$]")
	# ax.set_ylabel(r"$\varkappa_\nu, \sigma_\nu$ [cm$^2$ g$^{-1}$]")
	ax.set_xlim(9e-6,1e0)
	ax.set_ylim(1e-2,1e5)
	ax.tick_params(direction='in', which='both', pad=10, width=2, length=5)
	ax.tick_params(which='major', length=10)
	for side in ax.spines.values():
		side.set_linewidth(2)
	ax.loglog()
	ax.grid(linestyle=':')

	texts = ax.legend(loc="upper right").get_texts() 

	return texts


def add_text(ax):
	bbox = dict(boxstyle="round", fc="white", ec="gray", alpha=0.6)
	txt = r"$-\ \varkappa_\nu^{\rm abs}$"+"\n"+r"-- $\varkappa_\nu^{\rm sca}$"
	# txt = r"$-\ \varkappa_\nu$"+"\n"+r"-- $\sigma_\nu$"
	ax.text(0.03, 0.03, txt, va="bottom", ha="left", transform=ax.transAxes, alpha=0, bbox=bbox, rasterized=True)
	ax.text(0.03, 0.03, txt, va="bottom", ha="left", transform=ax.transAxes)

def add_legend2(ax):
	ax2 = ax.twinx()
	ax2.axis("off")

	ax2.plot(-1, -1, "-",  color="black", label=r"$\varkappa_\nu^{\rm abs}$")
	ax2.plot(-1, -1, "--", color="black", label=r"$\varkappa_\nu^{\rm sca}$")

	ax2.legend(loc="lower left")


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
name = "kappaNu_from_lambda"
run_comm(f"mkdir {name}", pwd_path)
res = "./" + name


plt.rcParams.update({'font.size': 30})


for sca in ["yes"]:

	pdf = PdfPages(f"{res}/{name}_amax_sca{sca}.pdf")

	plt.figure(figsize=(12,10))
	ax = plt.subplot(111)

	for i, amax_micron in enumerate(np.sort(np.append(amax_arr_micron, amax0_micron))):
		meta_name = f"p{p0}f{frac10}amax{amax_micron}sca{sca}"

		# if amax_micron == amax0_micron: k = i

		lambd, kappa, sigma = get_data(meta_name)

		if amax_micron == amax0_micron:
			ax.plot(lambd, kappa, "-", color=f"C{i}", label=rf"$\mathbf{{a_{{max}} = {amax_micron}\ \mu m}}$")
		else: 
			ax.plot(lambd, kappa, "-", color=f"C{i}", label=rf"$a_{{\rm max}} = {amax_micron}\ \mu \rm m$")
		# ax.plot(lambd, kappa, "-", color=f"C{i}",  label=r"a$_{\rm max} =$"+f"{amax_micron} $\mu$m")
		ax.plot(lambd, sigma, "--", color=f"C{i}")

	texts = set_ax(ax)
	# texts[k].set_weight("bold")

	plt.title(rf"$p = {p0},\ f_{{\rm Si}} = {frac10}$")
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

	for i, frac1 in enumerate(np.sort(np.append(frac1_arr, frac10))):
		meta_name = f"p{p0}f{frac1}amax{amax0_micron}sca{sca}"

		# if frac1 == frac10: k = i

		lambd, kappa, sigma = get_data(meta_name)

		if frac1 == frac10:
			ax.plot(lambd, kappa, "-", color=f"C{i}", label=rf"$\mathbf{{f_{{Si}} = {frac1}}}$")
		else:
			ax.plot(lambd, kappa, "-", color=f"C{i}", label=rf"$f_{{\rm Si}} = {frac1}$")
		# ax.plot(lambd, kappa, "-", color=f"C{i}", label=r"frac$_{\rm Si} =$"+f"{frac1}")
		ax.plot(lambd, sigma, "--", color=f"C{i}")

	texts = set_ax(ax)
	# texts[k].set_weight("bold")

	plt.title(rf"$p = {p0},\ a_{{\rm max}} = {amax0_micron}\ \mu \rm m$")
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

	for i, p in enumerate(np.sort(np.append(p_arr, p0))):
		meta_name = f"p{p}f{frac10}amax{amax0_micron}sca{sca}"

		# if p == p0: k = i

		lambd, kappa, sigma = get_data(meta_name)

		if p == p0:
			ax.plot(lambd, kappa, "-", color=f"C{i}", label=rf"$\mathbf{{p = {p}}}$")
		else:
			ax.plot(lambd, kappa, "-", color=f"C{i}", label=rf"$p = {p}$")
		# ax.plot(lambd, kappa, "-", color=f"C{i}", label=f"p = {p}")
		ax.plot(lambd, sigma, "--", color=f"C{i}")

	texts = set_ax(ax)
	# texts[k].set_weight("bold")

	# plt.title(f"fracSi = {frac10}, amax = {amax0_micron} $\mu$m")
	plt.title(rf"$f_{{\rm Si}} = {frac10},\ a_{{\rm max}} = {amax0_micron}\ \mu \rm m$")
	# add_text(ax)
	add_legend2(ax)

	plt.tight_layout()
	pdf.savefig()
	plt.close()
	
	pdf.close()

