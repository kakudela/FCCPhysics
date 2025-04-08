intLumi        = 1.0 # assume histograms are scaled in previous step
outputDir      = "output/ZmumuH/combine/"
mc_stats       = False
rebin          = 1

# get histograms from histmaker step
inputDir       = "output/ZmumuH/histmaker/"
# You can also use the pre-generated files
# /ceph/submit/data/group/fcc/ee/tutorials/FNAL/ZmumuH/histmaker/ MIT
# /ceph/submit/data/group/fcc/ee/tutorials/FNAL/ZmumuH/histmaker/ CERN


sig_procs = {'sig':['wzp6_ee_mumuH_Hbb_ecm240', 'wzp6_ee_mumuH_Hcc_ecm240', 'wzp6_ee_mumuH_Hss_ecm240', 'wzp6_ee_mumuH_Hgg_ecm240', 'wzp6_ee_mumuH_Haa_ecm240', 'wzp6_ee_mumuH_HZa_ecm240', 'wzp6_ee_mumuH_HWW_ecm240', 'wzp6_ee_mumuH_HZZ_ecm240', 'wzp6_ee_mumuH_Hmumu_ecm240', 'wzp6_ee_mumuH_Htautau_ecm240']}

# for the backgrounds, we merge all processes
bkg_procs = {'bkg':['wzp6_ee_mumu_ecm240', 'wzp6_ee_tautau_ecm240', 'p8_ee_WW_mumu_ecm240', 'p8_ee_ZZ_ecm240']}



categories = ["zmumu"]
hist_names = ["zmumu_recoil_m_final"]


systs = {}

# 2 % uncertainty on background normalization
#systs['bkg_norm'] = {
#    'type': 'lnN',
#    'value': 1.02,
#    'procs': ['bkg'],
#}

