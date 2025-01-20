import ROOT

# global parameters
intLumi        = 1. # assume histograms are scaled in previous step
intLumiLabel   = "L = 44.84 pb^{-1}"
ana_tex        = 'e^{+}e^{-} #rightarrow Z #rightarrow #mu^{#plus}#mu^{#minus}'
energy         = 91.2
collider       = 'FCC-ee'
formats        = ['png','pdf']

outdir         = '/home/submit/jaeyserm/public_html/fccee/tutorials/z_mumu_xsec/'
inputDir       = './output/tutorials/z_mumu_xsec/' 

plotStatUnc    = True



procs = {}
procs['signal'] = {'mumu':['wzp6_ee_mumu_ecm91p2']}
procs['backgrounds'] =  {
    'gaga':['p8_ee_gaga_mumu_ecm91p2'],
    'tautau':['wzp6_ee_tautau_ecm91p2'],
    }


colors = {}
colors['mumu'] = ROOT.kRed
colors['gaga'] = ROOT.kBlue+1
colors['tautau'] = ROOT.kRed+2

legend = {}
legend['mumu'] = "#mu^{#plus}#mu^{#minus}"
legend['gaga'] = "e^{#plus}e^{#minus}qq"
legend['tautau'] = "#tau^{#plus}#tau^{#minus}"

hists = {}


hists["muons_all_costheta"] = {
    "output":   "muons_all_costheta",
    "logy":     True,
    "stack":    True,
    "rebin":    1,
    "xmin":     -1,
    "xmax":     1,
    "ymin":     0.1,
    "ymax":     1e8,
    "xtitle":   "cos(#theta)",
    "ytitle":   "Events",
}

hists["muon_max_p_norm"] = {
    "output":   "muon_max_p_norm",
    "logy":     True,
    "stack":    True,
    "rebin":    1,
    "xmin":     0,
    "xmax":     2,
    "ymin":     0.1,
    "ymax":     1e8,
    "xtitle":   "p(#mu_{max})/E_{beam}",
    "ytitle":   "Events",
}

hists["acolinearity"] = {
    "output":   "acolinearity",
    "logy":     True,
    "stack":    True,
    "rebin":    1,
    "xmin":     0,
    "xmax":     1,
    "ymin":     1e-2,
    "ymax":     1e8,
    "xtitle":   "Acolinearity (rad)",
    "ytitle":   "Events",
}


hists["invariant_mass"] = {
    "output":   "invariant_mass",
    "logy":     True,
    "stack":    True,
    "rebin":    10,
    "xmin":     0,
    "xmax":     150,
    "ymin":     1e-3,
    "ymax":     1e6,
    "xtitle":   "Invariant mass (GeV)",
    "ytitle":   "Events",
}
