

intLumi        = 1.0 # assume histograms are scaled in previous step
outputDir      = "./output/tutorials/zh_xsec/combine"
mc_stats       = True
rebin          = 1

# get histograms from histmaker step
inputDir       = "./output/tutorials/zh_xsec/"


sig_procs = {'sig':['wzp6_ee_mumuH_ecm240']}
bkg_procs = {'bkg':['p8_ee_WW_mumu_ecm240', 'p8_ee_ZZ_ecm240', 'wzp6_ee_mumu_ecm240']}


categories = ["zmumu"]
hist_names = ["zmumu_recoil_m_final"]


systs = {}

# 10 % uncertainty on background normalization
systs['bkg_norm'] = {
    'type': 'lnN',
    'value': 1.10,
    'procs': ['bkg'],
}
