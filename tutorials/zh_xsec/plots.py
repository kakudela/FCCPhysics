import ROOT

# global parameters
intLumi        = 1. # assume histograms are scaled in previous step
intLumiLabel   = "L = 10.8 ab^{-1}"
ana_tex        = 'e^{+}e^{-} #rightarrow ZH #rightarrow #mu^{+}#mu^{-} + X'
energy         = 240.0
collider       = 'FCC-ee'
formats        = ['png','pdf']

outdir         = '/home/submit/jaeyserm/public_html/fccee/tutorials/zh_xsec/'
inputDir       = './output/tutorials/zh_xsec/' 

plotStatUnc    = True


procs = {}
procs['signal'] = {'ZH':['wzp6_ee_mumuH_ecm240']}
procs['backgrounds'] =  {
    'WW':['p8_ee_WW_mumu_ecm240'],
    'ZZ':['p8_ee_ZZ_ecm240'],
    'Zgamma':['wzp6_ee_mumu_ecm240']
    }

colors = {}
colors['ZH'] = ROOT.kRed
colors['WW'] = ROOT.kBlue+1
colors['ZZ'] = ROOT.kGreen+2
colors['Zgamma'] = ROOT.kRed+2

legend = {}
legend['ZH'] = 'ZH'
legend['WW'] = 'WW'
legend['ZZ'] = 'ZZ'
legend['Zgamma'] = 'Z/#gamma*'


hists = {}


hists["zmumu_recoil_m_final"] = {
    "output":   "zmumu_recoil_m_final",
    "logy":     False,
    "stack":    True,
    "rebin":    1,
    "xmin":     120,
    "xmax":     140,
    "ymin":     0,
    "ymax":     2500,
    "xtitle":   "Recoil (GeV)",
    "ytitle":   "Events / 100 MeV",
}

hists["zmumu_recoil_m"] = {
    "output":   "zmumu_recoil_m",
    "logy":     False,
    "stack":    False,
    "rebin":    1,
    "xmin":     100,
    "xmax":     150,
    "ymin":     0,
    "ymax":     15000,
    "xtitle":   "Recoil (GeV)",
    "ytitle":   "Events / 100 MeV",
}

hists["zmumu_p"] = {
    "output":   "zmumu_p",
    "logy":     False,
    "stack":    False,
    "rebin":    1,
    "xmin":     0,
    "xmax":     120,
    "ymin":     0,
    "ymax":     100000,
    "xtitle":   "p(#mu^{#plus}#mu^{#minus}) (GeV)",
    "ytitle":   "Events ",
    "scaleSig": 10
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
