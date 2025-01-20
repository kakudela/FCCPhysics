import ROOT

# global parameters
intLumi        = 1. # assume histograms are scaled in previous step
intLumiLabel   = "L = 10.8 ab^{-1}"
ana_tex        = 'e^{+}e^{-} #rightarrow ZH #rightarrow #mu^{+}#mu^{-} + X'
energy         = 240.0
collider       = 'FCC-ee'
formats        = ['png','pdf']

outdir         = '/home/submit/jaeyserm/public_html/fccee/h_aa/plots/'
inputDir       = './output/h_aa/histmaker/' 

plotStatUnc    = False
leg_position = [0.8, 0.7, 0.96, 0.88]

procs = ["", "", ""]

procs = {}
procs['signal'] = {'ZH':['wzp6_ee_nunuH_Haa_ecm240', 'wzp6_ee_eeH_Haa_ecm240', 'wzp6_ee_tautauH_Haa_ecm240', 'wzp6_ee_ccH_Haa_ecm240', 'wzp6_ee_bbH_Haa_ecm240', 'wzp6_ee_qqH_Haa_ecm240', 'wzp6_ee_ssH_Haa_ecm240', 'wzp6_ee_mumuH_Haa_ecm240']}
procs['backgrounds'] =  {
    'Zgamma':['kkmcee_ee_uu_ecm240', 'kkmcee_ee_dd_ecm240', 'kkmcee_ee_cc_ecm240', 'kkmcee_ee_ss_ecm240', 'kkmcee_ee_bb_ecm240', 'kkmcee_ee_tautau_ecm240', 'kkmcee_ee_mumu_ecm240'],
    'gaga':['wzp6_ee_gammagamma_ecm240'],
    }


colors = {}
colors['ZH'] = ROOT.kRed
colors['WW'] = ROOT.kBlue+1
colors['gaga'] = ROOT.kGreen+2
colors['Zgamma'] = ROOT.kRed+2

legend = {}
legend['ZH'] = 'ZH'
legend['WW'] = 'WW'
legend['gaga'] = '#gamma#gamma'
legend['Zgamma'] = 'Z/#gamma* #rightarrow f#bar{f}#gamma(#gamma)'


hists = {}

hists["photons_all_p"] = {
    "output":   "photons_all_p",
    "logy":     True,
    "stack":    False,
    "rebin":    1,
    "xmin":     0,
    "xmax":     150,
    "ymin":     1e-1,
    "ymax":     -1,
    "xtitle":   "Photons momentum (GeV)",
    "ytitle":   "Events",
}

hists["photons_all_costheta"] = {
    "output":   "photons_all_costheta",
    "logy":     True,
    "stack":    False,
    "rebin":    1,
    "xmin":     -1,
    "xmax":     1,
    "ymin":     1e-1,
    "ymax":     -1,
    "xtitle":   "#cos#theta photons",
    "ytitle":   "Events",
}

hists["photons_all_theta"] = {
    "output":   "photons_all_theta",
    "logy":     True,
    "stack":    False,
    "rebin":    1,
    "xmin":     0,
    "xmax":     3.15,
    "ymin":     1e-1,
    "ymax":     -1,
    "xtitle":   "#theta photons (rad0",
    "ytitle":   "Events",
    "scaleSig": 1
}


hists["photons_sel_p_costheta"] = {
    "output":   "photons_sel_p_costheta",
    "logy":     True,
    "stack":    False,
    "rebin":    1,
    "xmin":     -1,
    "xmax":     1,
    "ymin":     1e-1,
    "ymax":     -1,
    "xtitle":   "#cos#theta photons",
    "ytitle":   "Events",
}

hists["photons_sel_p_theta"] = {
    "output":   "photons_sel_p_theta",
    "logy":     True,
    "stack":    False,
    "rebin":    1,
    "xmin":     0,
    "xmax":     3.15,
    "ymin":     1e-1,
    "ymax":     -1,
    "xtitle":   "#theta photons (rad0",
    "ytitle":   "Events",
    "scaleSig": 1
}




hists["zmumu_m"] = {
    "output":   "zmumu_m",
    "logy":     False,
    "stack":    True,
    "rebin":    1,
    "xmin":     0,
    "xmax":     250,
    "ymin":     0,
    "ymax":     3000,
    "xtitle":   "m(#mu^{#plus}#mu^{#minus}) (GeV)",
    "ytitle":   "Events ",
}

hists["cosThetaMiss_cut4"] = {
    "output":   "cosThetaMiss_cut4",
    "logy":     True,
    "stack":    False,
    "rebin":    10,
    "xmin":     0,
    "xmax":     1,
    "ymin":     10,
    "ymax":     100000,
    "xtitle":   "cos(#theta_{miss})",
    "ytitle":   "Events ",
}


hists["cutFlow"] = {
    "output":   "cutFlow",
    "logy":     True,
    "stack":    False,
    "xmin":     0,
    "xmax":     6,
    "ymin":     1e4,
    "ymax":     1e10,
    "xtitle":   ["All events", "#geq 1 #mu^{#pm} + ISO", "#geq 2 #mu^{#pm} + OS", "86 < m_{#mu^{+}#mu^{#minus}} < 96", "20 < p_{#mu^{+}#mu^{#minus}} < 70", "|cos#theta_{miss}| < 0.98", "120 < m_{rec} < 140"],
    "ytitle":   "Events ",
    "scaleSig": 10
}
