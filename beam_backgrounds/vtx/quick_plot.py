import os, glob, ROOT, itertools
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

datasets = [
    ("/ceph/submit/data/user/k/kudela/beam_backgrounds/altered_gp_pairs/"
     "FCCee_Z_4IP_04may23_FCCee_Z128_edge",
     "128_original_mapped"),

    ("/ceph/submit/data/user/k/kudela/beam_backgrounds/altered_gp_pairs/"
     "FCCee_Z_4IP_04may23_FCCee_Z128_grids1_testGrids_edge",
     "128_testGrids_fix_mapped"),
]

outDir_base = "/home/submit/kudela/fccee/FCCPhysics/beam_backgrounds/vtx/quick_mc_plots"
compareDir  = os.path.join(outDir_base, "compare")
os.makedirs(compareDir, exist_ok=True)

HX, HY, HZ = 1.32555, 0.00186, 25.4
buf = 1.1
x_lo, x_hi = -buf*HX, buf*HX
y_lo, y_hi = -buf*HY, buf*HY
z_lo, z_hi = -buf*HZ, buf*HZ
nb  = 400

h1d  = {lab: {} for _, lab in datasets}
h2d  = {lab: {} for _, lab in datasets}
stats= {lab: dict(sumE=0,sumPx=0,sumPy=0,sumPz=0,n=0) for _,lab in datasets}

def new_h1(name, title, bins, lo, hi):
    return ROOT.TH1F(name, title, bins, lo, hi)

for inDir, lab in datasets:

    outDir = os.path.join(outDir_base, lab)
    os.makedirs(outDir, exist_ok=True)

    # 1-D hists
    hx  = new_h1(f"h_x_{lab}",  ";x (mm);Entries", nb, x_lo, x_hi)
    hy  = new_h1(f"h_y_{lab}",  ";y (mm);Entries", nb, y_lo, y_hi)
    hz  = new_h1(f"h_z_{lab}",  ";z (mm);Entries", nb, z_lo, z_hi)
    hE  = new_h1(f"h_E_{lab}",  ";|E| (GeV);Entries", 400, 0, 50)
    hE_MeV = new_h1(f"h_E_MeV_{lab}", ";|E| (MeV);Entries", 400, 0, 10)
    hpx = new_h1(f"h_px_{lab}",";|beta_x|;Entries", 400, 0, 1.1)
    hpy = new_h1(f"h_py_{lab}",";|beta_y|;Entries", 400, 0, 1.1)
    hpz = new_h1(f"h_pz_{lab}",";|beta_z|;Entries", 400, 0, 1.1)
    hb  = new_h1(f"h_beta_{lab}",";|beta|;Entries", 400, 0, 1.5)

    # 2-D hists
    h_xy = ROOT.TH2F(f"h_xy_{lab}",";x (mm);y (mm)", nb, x_lo, x_hi, nb, y_lo, y_hi)
    h_xz = ROOT.TH2F(f"h_xz_{lab}",";x (mm);z (mm)", nb, x_lo, x_hi, nb, z_lo, z_hi)
    h_yz = ROOT.TH2F(f"h_yz_{lab}",";y (mm);z (mm)", nb, y_lo, y_hi, nb, z_lo, z_hi)

    h1d[lab] = dict(E=hE, EMeV=hE_MeV, px=hpx, py=hpy, pz=hpz, beta=hb,
                    x=hx, y=hy, z=hz)
    h2d[lab] = dict(xy=h_xy, xz=h_xz, yz=h_yz)

    st = stats[lab]

    for fp in glob.glob(os.path.join(inDir, "*.pairs")):
        with open(fp) as fh:
            for ln in fh:
                if not ln.strip() or ln.startswith("E"):
                    continue
                try:
                    E, px, py, pz, x_nm, y_nm, z_nm, *_ = map(float, ln.split())
                    x_mm, y_mm, z_mm = x_nm*1e-6, y_nm*1e-6, z_nm*1e-6

                    hx.Fill(x_mm); hy.Fill(y_mm); hz.Fill(z_mm)
                    h_xy.Fill(x_mm, y_mm); h_xz.Fill(x_mm, z_mm); h_yz.Fill(y_mm, z_mm)

                    absE = abs(E)
                    hE.Fill(absE); hE_MeV.Fill(absE*1e3)
                    hpx.Fill(abs(px)); hpy.Fill(abs(py)); hpz.Fill(abs(pz))
                    hb.Fill((px*px + py*py + pz*pz)**0.5)

                    st["sumE"]  += absE
                    st["sumPx"] += abs(px)
                    st["sumPy"] += abs(py)
                    st["sumPz"] += abs(pz)
                    st["n"]     += 1
                except ValueError:
                    pass

    nFiles = len(glob.glob(os.path.join(inDir, "*.pairs")))
    if nFiles > 0:
        scale = 1.0 / nFiles
        for h in h1d[lab].values(): h.Scale(scale)
        for h in h2d[lab].values(): h.Scale(scale)

    # save per-dataset 2-D maps
    cv = ROOT.TCanvas("", "", 1000, 1000)
    cv.SetLeftMargin(0.12); cv.SetBottomMargin(0.11)
    for key, h in h2d[lab].items():
        h.Draw("COLZ")
        cv.SaveAs(os.path.join(outDir, f"{key}.png"))
        cv.Clear()

color_cycle = [ROOT.kBlue+1, ROOT.kRed+1, ROOT.kGreen+2,
               ROOT.kMagenta+1, ROOT.kCyan+1, ROOT.kOrange+7]

variables = [
    ("x",   ";x (mm);Entries"),
    ("y",   ";y (mm);Entries"),
    ("z",   ";z (mm);Entries"),
    ("E",   ";|E| (GeV);Entries"),
    ("EMeV",";|E| (MeV);Entries"),
    ("px",  ";|beta_x|;Entries"),
    ("py",  ";|beta_y|;Entries"),
    ("pz",  ";|beta_z|;Entries"),
    ("beta",";|beta|;Entries"),
]

cv = ROOT.TCanvas("", "", 1000, 1000)
cv.SetLeftMargin(0.12); cv.SetBottomMargin(0.11)

for var, axisTitle in variables:
    ymax = max(h1d[lab][var].GetMaximum() for _, lab in datasets)
    leg = ROOT.TLegend(0.58,0.75,0.88,0.88)
    leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.02)
    first = True
    for col,( _, lab ) in zip(itertools.cycle(color_cycle), datasets):
        h = h1d[lab][var]
        h.SetLineColor(col); h.SetLineWidth(2)
        h.GetXaxis().SetTitle(axisTitle.split(";")[1])
        h.GetYaxis().SetTitle("entries / file")
        h.SetMaximum(ymax*1.10)
        h.Draw("HIST" if first else "HIST SAME")
        first = False
        leg.AddEntry(h, lab, "l")
    leg.Draw()
    cv.SaveAs(os.path.join(compareDir, f"{var}.png"))
    cv.Clear()

for _, lab in datasets:
    st = stats[lab]
    if st["n"]:
        print(f"\nDataset: {lab}")
        print(f"  <|E|>      = {st['sumE']/st['n']:.3f} GeV")
        print(f"  <|beta_x|> = {st['sumPx']/st['n']:.4f}")
        print(f"  <|beta_y|> = {st['sumPy']/st['n']:.4f}")
        print(f"  <|beta_z|> = {st['sumPz']/st['n']:.4f}")

for _, lab in datasets:
    outDir = os.path.join(outDir_base, lab)
    fout   = ROOT.TFile(os.path.join(outDir, "output.root"), "RECREATE")
    for h in h1d[lab].values(): h.Write()
    for h in h2d[lab].values(): h.Write()
    fout.Close()
