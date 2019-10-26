#!/usr/bin/env python

"""
Bayu Wilson, 2018-2019, University of Washington Astronomy

Usage

Purpose

"""
import os

os.system("./main.py")
os.chdir("../plot/")
os.system("./plot_mf.py")
os.system("./plot_pk.py")
os.chdir("../pipeline")

# os.system("./main.py")
# os.chdir("../plot/")
# # os.system("./plot_mf.py")
# os.system("./plot_specific.py")
# os.chdir("../pipeline")


#os.system('./bootstrap_pk.py')
# os.chdir("../plot/")
# os.system("./plot_ratio_pk.py")
# os.system("./plot_nocorr_ratio.py")
# os.system("./plot_bootstrap.py")
# os.chdir("../pipeline")


# os.system("./findmf.py mocks lyb_nocorr 1")
# os.system("./findpk.py mocks lyb_nocorr 1")

# os.system("./findmf.py mocks lyb_nocorr_n5000 1")
# os.system("./findpk.py mocks lyb_nocorr_n5000 1")


# Mean Flux
### XQ-100 Observations
# os.system("./findmf.py obs uncorr 1")
# os.system("./findmf.py obs corrN 1")
# os.system("./findmf.py obs corrNR 1")
# os.system("./findmf.py obs corrR 1")

# os.system("./findmf.py mocks nocorr 1")
# os.system("./findmf.py mocks wNR 1")
# os.system("./findmf.py mocks wN 1")
# os.system("./findmf.py mocks wR 1")
# os.system("./findmf.py mocks wN_SN10  1")
# os.system("./findmf.py mocks wR_n5000 1")
# os.system("./findmf.py mocks wN_n5000 1")
# os.system("./findmf.py mocks wNR_n5000 1")
#
# # os.system("./findmf.py mocks lyb_wN 1")
# # os.system("./findmf.py mocks lyb_wR 1")
# # os.system("./findmf.py mocks lyb_wNR 1")
# os.system("./findmf.py mocks lyb_nocorr 1")
# os.system("./findmf.py mocks lyb_nocorr_n5000 1")
# os.system("./findmf.py mocks noB 1")
# os.system("./findmf.py mocks noB_onlyA 1")

# os.system("./findmf.py mocks lyb_nocorr 1")
# os.system("./findmf.py mocks lyb_nocorr_n5000 1")
# os.system("./findmf.py mocks lyb_wN_n5000 1")
# # os.system("./findmf.py mocks lyb_wR_n5000 1") #PROBLEM
# os.system("./findmf.py mocks lyb_wNR_n5000 1")
# os.system("./findmf.py mocks noB 1")
# os.system("./findmf.py mocks noB_onlyA 1")


# Power spectrum

### XQ-100 OBSERVATIONS
# os.system("./findpk.py obs uncorr 1")
# os.system("./findpk.py obs corrN 1")
# os.system("./findpk.py obs corrR 1")
# os.system("./findpk.py obs corrNR 1")

### LYB MOCKS n5000
# os.system("./findpk.py mocks lyb_nocorr_n5000 1")
# os.system("./findpk.py mocks lyb_wN_n5000 1")
# # os.system("./findpk.py mocks lyb_wR_n5000 1") #PROBLEM
# os.system("./findpk.py mocks lyb_wNR_n5000 1")



# os.system("./findpk.py obs uncorr 1")
# os.system("./findpk.py mocks nocorr 1")
# os.system("./findpk.py mocks wNR 1")
# os.system("./findpk.py mocks wN 1")
# os.system("./findpk.py mocks wR 1")
# os.system("./findpk.py mocks wN_SN10  1")
# os.system("./findpk.py mocks wR_n5000 1")
# os.system("./findpk.py mocks wN_n5000 1")
# os.system("./findpk.py mocks wNR_n5000 1")

# os.system("./findmf.py mocks lyb_wN 1")
# os.system("./findmf.py mocks lyb_wR 1")
# os.system("./findmf.py mocks lyb_wNR 1")
# os.system("./findpk.py mocks lyb_nocorr 1")
# os.system("./findpk.py mocks lyb_wN_n5000 1")
# os.system("./findpk.py mocks lyb_wR_n5000 1") #PROBLEM
# os.system("./findpk.py mocks lyb_wNR_n5000 1")
# os.system("./findpk.py mocks noB 1")
# os.system("./findpk.py mocks noB_onlyA 1")

# os.system("./findpk.py mocks lyb_nocorr 1")



#
# # Bootstrap
# # os.system("./bootstrap_lyb.py mocks lyb_nocorr 1 5")
# # os.system("./bootstrap_lyb.py obs uncorr 1 5")
# # os.system("./bootstrap_lyb.py obs corrN 1 5")
# # os.system("./bootstrap_lyb.py obs corrR 1 5")
# os.system("./bootstrap_lyb.py obs corrNR 1 1000")
#
#
# # Making plots
# os.chdir("plots")
# # os.system("./pk_tot_plot.py mocks lyb_nocorr 1")
# # os.system("./pk_tot_plot.py obs unscorr 1")
# # os.system("./pk_tot_plot.py obs corrN 1")
# # os.system("./pk_tot_plot.py obs corrR 1")
# os.system("./pk_tot_plot.py obs corrNR 1")
# os.chdir("../")
