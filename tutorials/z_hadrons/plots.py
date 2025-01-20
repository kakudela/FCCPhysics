import ROOT

# global parameters
intLumi        = 1. # assume histograms are scaled in previous step
intLumiLabel   = "L = 44.84 pb^{-1}"
ana_tex        = 'e^{+}e^{-} #rightarrow Z #rightarrow hadrons'
energy         = 91.2
collider       = 'FCC-ee'
formats        = ['png','pdf']

outdir         = '/home/submit/jaeyserm/public_html/fccee/tutorials/z_hadrons/'
inputDir       = './output/tutorials/z_hadrons/' 

plotStatUnc    = True



procs = {}
procs['signal'] = {'qq':['wz3p6_ee_qq_ecm91p2']}
procs['backgrounds'] =  {
    'gaga':['p8_ee_gaga_qq_ecm91p2'],
    'ee':['wzp6_ee_ee_Mee_5_150_ecm91p2'],
    'mumu':['wzp6_ee_mumu_ecm91p2'],
    'tautau':['wzp6_ee_tautau_ecm91p2'],
    }


colors = {}
colors['qq'] = ROOT.kRed
colors['gaga'] = ROOT.kBlue+1
colors['ee'] = ROOT.kGreen+2
colors['mumu'] = ROOT.kRed+2
colors['tautau'] = ROOT.kRed

legend = {}
legend['qq'] = "qq"
legend['gaga'] = "e^{#plus}e^{#minus}qq"
legend['ee'] = "e^{#plus}e^{#minus}"
legend['mumu'] = "#mu^{#plus}#mu^{#minus}"
legend['tautau'] = "#tau^{#plus}#tau^{#minus}"

hists = {}


hists["visible_energy_norm_nOne"] = {
    "output":   "visible_energy_norm_nOne",
    "logy":     True,
    "stack":    True,
    "rebin":    1,
    "xmin":     0.2,
    "xmax":     2.4,
    "ymin":     0.1,
    "ymax":     1e8,
    "xtitle":   "E_{vis}/#sqrt{s}",
    "ytitle":   "Events",
}

hists["thrust_costheta_nCut"] = {
    "output":   "thrust_costheta_nCut",
    "logy":     True,
    "stack":    True,
    "rebin":    1,
    "xmin":     0,
    "xmax":     1,
    "ymin":     10,
    "ymax":     1e7,
    "xtitle":   "|cos(#theta_{t})|",
    "ytitle":   "Events",
}

hists["thrust_magn_nCut"] = {
    "output":   "thrust_magn_nCut",
    "logy":     True,
    "stack":    True,
    "rebin":    1,
    "xmin":     0,
    "xmax":     1,
    "ymin":     10,
    "ymax":     1e7,
    "xtitle":   "Thrust magnitude",
    "ytitle":   "Events",
}


hists["energy_imbalance_long_nOne"] = {
    "output":   "energy_imbalance_long_nOne",
    "logy":     True,
    "stack":    True,
    "rebin":    1,
    "xmin":     0,
    "xmax":     1,
    "ymin":     0.1,
    "ymax":     1e7,
    "xtitle":   "E_{long}/E_{vis}",
    "ytitle":   "Events",
}

hists["energy_imbalance_trans_nOne"] = {
    "output":   "energy_imbalance_trans_nOne",
    "logy":     True,
    "stack":    True,
    "rebin":    1,
    "xmin":     0,
    "xmax":     1,
    "ymin":     0.1,
    "ymax":     1e7,
    "xtitle":   "E_{trans}/E_{vis}",
    "ytitle":   "Events",
}

hists["RP_no_barrel_nOne"] = {
    "output":   "RP_no_barrel_nOne",
    "logy":     True,
    "stack":    True,
    "rebin":    1,
    "xmin":     0,
    "xmax":     120,
    "ymin":     0.1,
    "ymax":     1e7,
    "xtitle":   "Number of particles (barrel)",
    "ytitle":   "Events",
}


hists["RP_no_endcap_nOne"] = {
    "output":   "RP_no_endcap_nOne",
    "logy":     True,
    "stack":    True,
    "rebin":    1,
    "xmin":     0,
    "xmax":     120,
    "ymin":     0.1,
    "ymax":     1e7,
    "xtitle":   "Number of particles (endcap)",
    "ytitle":   "Events",
}


