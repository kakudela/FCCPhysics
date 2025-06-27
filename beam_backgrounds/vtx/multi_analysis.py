#!/usr/bin/env python3
# multi_analysis.py  – build graphs & write one ROOT file / dataset
#   * all histograms are **raw counts** (normalisation happens in the plotting scripts)
#   * threads auto-detected or set with --nThreads
#   * edit the datasets list at the end to add / remove samples
# ------------------------------------------------------------------------------

import ROOT, argparse, glob, os, time, logging
ROOT.gInterpreter.Declare('#include "functions.h"')

# ------- CLI -------
ap = argparse.ArgumentParser()
ap.add_argument('--nThreads',  type=int, default=None)
ap.add_argument('--maxFiles',  type=int, default=-1)
ap.add_argument('--outputDir', type=str,  default='output')
args = ap.parse_args()

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
log = logging.getLogger('multi')

# ------- detector geometry & binning -------
layer_radii = [14., 36., 58.]
max_z   = 110
bins_z  = (max_z // 2, -max_z, max_z)       # 55 bins, −110 → 110 mm
bins_th = (36, 0, 180)                      # 5-deg bins
bins_ph = (72, -180, 180)                   # 5-deg bins
bins_layer = (len(layer_radii), 0, len(layer_radii))

# ------------------------------------------------------------------------------
def build_graph(df):
    radii_cpp = '{' + ','.join(f'{r}f' for r in layer_radii) + '}'

    # ---- columns ----
    df = (
        df.Define('h_r',      'getSimHitPosition_r(VertexBarrelCollection)')
          .Define('h_z',      'getSimHitPosition_z(VertexBarrelCollection)')
          .Define('h_th',     'getSimHitPosition_theta(VertexBarrelCollection,true)')
          .Define('h_ph',     'getSimHitPosition_phi(VertexBarrelCollection,true)')
          .Define('h_E',      'getSimHitEDep(VertexBarrelCollection)')
          .Define('h_len',    'getSimHitPathLength(VertexBarrelCollection)')
          .Define('h_dEdx',   'h_E/h_len')
          .Define('h_cell',   'getSimHitCellID(VertexBarrelCollection)')
          .Define('h_layer',  f'getSimHitLayer(h_r,{radii_cpp})')
          .Define('h_m',      'filterSimHitsPrimary(VertexBarrelCollection)')
          .Define('h_L0',     'h_layer==0 && h_m')
          .Define('h_z_L0',   'h_z[h_L0]')
          .Define('h_th_L0',  'h_th[h_L0]')
          .Define('h_ph_L0',  'h_ph[h_L0]')
          .Define('h_cell_L0','h_cell[h_L0]')
          .Define('h_layer_m','h_layer[h_m]')
          .Define('h_dEdx_noSec','h_dEdx[h_m]')                 # primaries only
        # MC truth
          .Define('mc_z',     'getMCParticleVertex_z(MCParticles)')
          .Define('mc_r',     'getMCParticleVertex_r(MCParticles)')
          .Define('mc_th',    'getMCParticleVertex_theta(MCParticles,true)')
          .Define('mc_ph',    'getMCParticleVertex_phi(MCParticles,true)')
          .Define('mc_E',     'getMCParticleEnergy(MCParticles)')
          .Define('mc_m',     'filterGenStatus1(MCParticles)')
          .Define('mc_z_sel','mc_z[mc_m]').Define('mc_r_sel','mc_r[mc_m]')
          .Define('mc_th_sel','mc_th[mc_m]').Define('mc_ph_sel','mc_ph[mc_m]')
          .Define('mc_E_sel','mc_E[mc_m]')
          # momentum–theta (deg)
          .Define('mc_thetaMom',
                  'Vec_d o(mc_E.size());'                  # reuse size
                  'for(size_t i=0;i<o.size();++i){'
                  ' double px=MCParticles[i].momentum[0];'
                  ' double py=MCParticles[i].momentum[1];'
                  ' double pz=MCParticles[i].momentum[2];'
                  ' double p=std::sqrt(px*px+py*py+pz*pz);'
                  ' o[i]=p>0?std::acos(pz/p)*180.0/M_PI:0.;}'
                  'return o;')
    )

    # ---- histograms ----
    H=[]; add=H.append
    for n,v,b in [
        ('z','h_z',bins_z), ('theta','h_th',bins_th), ('phi','h_ph',bins_ph),
        ('z_L0','h_z_L0',bins_z), ('theta_L0','h_th_L0',bins_th),
        ('phi_L0','h_ph_L0',bins_ph), ('layer','h_layer',bins_layer)
    ]:
        add(df.Histo1D((n,'',*b),v))

    add(df.Histo2D(('z_phi','',*(bins_z+bins_ph)),'h_z','h_ph'))
    add(df.Histo2D(('theta_phi','',*(bins_th+bins_ph)),'h_th','h_ph'))
    add(df.Histo2D(('z_phi_L0','',*(bins_z+bins_ph)),'h_z_L0','h_ph_L0'))
    add(df.Histo2D(('theta_phi_L0','',*(bins_th+bins_ph)),'h_th_L0','h_ph_L0'))

    nL=len(layer_radii)
    add(df.Histo2D(('energy','',nL,0,nL,200,0,200),'h_layer','h_E'))
    add(df.Histo2D(('dEdx','',nL,0,nL,100,0,1000),'h_layer','h_dEdx'))
    add(df.Histo2D(('dEdx_noSec','',nL,0,nL,100,0,1000),'h_layer_m','h_dEdx_noSec'))

    add(df.Histo2D(('mc_z_r','',110,-110,110,200,0,110),'mc_z_sel','mc_r_sel'))
    add(df.Histo2D(('mc_z_phi','',*(bins_z+bins_ph)),'mc_z_sel','mc_ph_sel'))
    add(df.Histo2D(('mc_theta_phi','',*(bins_th+bins_ph)),'mc_th_sel','mc_ph_sel'))
    add(df.Histo1D(('mc_E','',100,0,100),'mc_E_sel'))
    add(df.Histo1D(('mc_thetaMom','',36,0,180),'mc_thetaMom'))

    occ = df.Histo1D(('cell_counts','',100000,0,100000),'h_cell_L0')
    nEv = df.Count()
    return H, occ, nEv, [nEv, occ] + H

# ------------------------------------------------------------------------------
def build_and_run(dsets):
    ROOT.EnableImplicitMT(args.nThreads or 0)
    log.info(f'Threads: {ROOT.GetThreadPoolSize()}')
    t0=time.time()

    graphs, triggers = [], []
    for name, path in dsets:
        files = glob.glob(f'{path}/*.root')
        if args.maxFiles>0: files = files[:args.maxFiles]
        if not files:
            log.warning(f'{name}: no files'); continue
        log.info(f'{name}: {len(files)} files')
        h,occ,nEv,tr = build_graph(ROOT.RDataFrame('events',files))
        graphs.append((name,h,occ,nEv)); triggers += tr

    log.info(f'Graphs built in {time.time()-t0:.1f} s')
    t_loop=time.time(); ROOT.RDF.RunGraphs(triggers)
    log.info(f'Event loop finished in {time.time()-t_loop:.1f} s')

    os.makedirs(args.outputDir,exist_ok=True)
    for name,h,occ,nEv in graphs:
        f=ROOT.TFile(f'{args.outputDir}/{name}.root','RECREATE')
        for hh in h: hh.GetValue().Write()
        occ.GetValue().Write()
        ROOT.TParameter(int)('nEvents',int(nEv.GetValue())).Write()
        f.Close()

    log.info(f'Total wall-time: {time.time()-t0:.1f} s')

# ------------------------------------------------------------------------------
def main():
    datasets=[
        ('256','/ceph/submit/data/user/k/kudela/beam_backgrounds/CLD_o2_v05/FCCee_Z_4IP_04may23_FCCee_Z256'),
        ('128','/ceph/submit/data/user/k/kudela/beam_backgrounds/CLD_o2_v05/FCCee_Z_4IP_04may23_FCCee_Z128')
    ]
    build_and_run(datasets)

if __name__=='__main__':
    main()





# import os, glob, time, pathlib, argparse, logging, math
# import ROOT

# parser = argparse.ArgumentParser()
# parser.add_argument("--nThreads",  type=int,  default=None)
# parser.add_argument("--maxFiles",  type=int,  default=-1)
# parser.add_argument("--outputDir", type=str,  default="output")

# parser.add_argument("--calculate", action="store_true")
# parser.add_argument("--plots",     action="store_true")
# parser.add_argument("--stats",     action="store_true")
# parser.add_argument("--log",       action="store_true")
# args = parser.parse_args()


# logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
# log = logging.getLogger("multi-analysis")


# if args.nThreads is not None:
#     ROOT.DisableImplicitMT()
#     ROOT.EnableImplicitMT(args.nThreads)
# else:
#     ROOT.EnableImplicitMT()
# log.info("Running with %d threads", ROOT.GetThreadPoolSize())


# ROOT.gROOT.SetBatch(True)
# ROOT.gStyle.SetOptStat(0)
# ROOT.gInterpreter.Declare('#include "functions.h"')


# def _canvas():
#     c = ROOT.TCanvas("c","c",800,700); c.SetRightMargin(0.15); return c
# def plot_hist(h,out,title):
#     c=_canvas(); h.SetTitle(title); h.Draw("HIST"); c.SaveAs(out); c.Close()
# def plot_2d(h,out,title,logz=False):
#     c=_canvas();
#     if logz:c.SetLogz(True)
#     h.SetTitle(title); h.Draw("COLZ"); c.SaveAs(out)
#     if logz:c.SetLogz(False); c.Close()


# # analysis constants (identical to analysis_2.py)
# layer_radii = [14., 36., 58.]
# max_z       = 110
# bins_theta  = (36, 0, 180)
# bins_phi    = (72, -180, 180)
# bins_z      = (max_z//2, -max_z, max_z)


# def my_graph(df, _):
#     # hit definitions
#     df = (df.Define("hits_r",     "getSimHitPosition_r(VertexBarrelCollection)")
#             .Define("hits_z",     "getSimHitPosition_z(VertexBarrelCollection)")
#             .Define("hits_theta", "getSimHitPosition_theta(VertexBarrelCollection,true)")
#             .Define("hits_phi",   "getSimHitPosition_phi(VertexBarrelCollection,true)")
#             .Define("hits_edep",  "getSimHitEDep(VertexBarrelCollection)")
#             .Define("hits_path",  "getSimHitPathLength(VertexBarrelCollection)")
#             .Define("hits_dedx",  "hits_edep/hits_path")
#             .Define("hits_cell",  "getSimHitCellID(VertexBarrelCollection)")
#             .Define("hits_layer", "getSimHitLayer(hits_r,{14.f,36.f,58.f})")
#             .Define("primary_mask","filterSimHitsPrimary(VertexBarrelCollection)")
#             .Define("mask_l0",    "hits_layer==0 && primary_mask")
#             .Define("hits_z_l0",     "hits_z[mask_l0]")
#             .Define("hits_theta_l0", "hits_theta[mask_l0]")
#             .Define("hits_phi_l0",   "hits_phi[mask_l0]")
#             .Define("hits_cell_l0",  "hits_cell[mask_l0]"))
#     # MC
#     df = (df.Define("mc_r",      "getMCParticleVertex_r(MCParticles)")
#             .Define("mc_z",      "getMCParticleVertex_z(MCParticles)")
#             .Define("mc_theta",  "getMCParticleVertex_theta(MCParticles,true)")
#             .Define("mc_phi",    "getMCParticleVertex_phi(MCParticles,true)")
#             .Define("mc_E",      "getMCParticleEnergy(MCParticles)")
#             .Define("mc_mask",   "filterGenStatus1(MCParticles)")
#             .Define("mc_z_sel",  "mc_z[mc_mask]")
#             .Define("mc_r_sel",  "mc_r[mc_mask]")
#             .Define("mc_theta_sel","mc_theta[mc_mask]")
#             .Define("mc_phi_sel","mc_phi[mc_mask]")
#             .Define("mc_E_sel",  "mc_E[mc_mask]"))

#     # histograms
#     H=[]
#     for n,v,b in [("z","hits_z",bins_z),
#                   ("theta","hits_theta",bins_theta),
#                   ("phi","hits_phi",bins_phi),
#                   ("z_layer0","hits_z_l0",bins_z),
#                   ("theta_layer0","hits_theta_l0",bins_theta),
#                   ("phi_layer0","hits_phi_l0",bins_phi)]:
#         H.append(df.Histo1D((n,"",*b),v))

#     H+= [df.Histo2D(("z_phi","",*(bins_z+bins_phi)),"hits_z","hits_phi"),
#          df.Histo2D(("theta_phi","",*(bins_theta+bins_phi)),"hits_theta","hits_phi"),
#          df.Histo2D(("z_phi_layer0","",*(bins_z+bins_phi)),"hits_z_l0","hits_phi_l0"),
#          df.Histo2D(("theta_phi_layer0","",*(bins_theta+bins_phi)),"hits_theta_l0","hits_phi_l0"),
#          df.Histo2D(("energy","",len(layer_radii),0,len(layer_radii),200,0,200),"hits_layer","hits_edep"),
#          df.Histo2D(("dedx","",len(layer_radii),0,len(layer_radii),100,0,1000),"hits_layer","hits_dedx"),
#          df.Histo2D(("mc_z_r","",110,-110,110,200,0,110),"mc_z_sel","mc_r_sel"),
#          df.Histo2D(("mc_theta_phi","",*(bins_theta+bins_phi)),"mc_theta_sel","mc_phi_sel"),
#          df.Histo1D(("mc_energy","",100,0,100),"mc_E_sel")]

#     occ = df.Histo1D(("cell_counts","",100_000,0,100_000),"hits_cell_l0")
#     nEv = df.Count()
#     return df, H, nEv, occ

# # core driver
# def build_and_run(datasets):
#     graphs = []
#     per_ds=[]
#     t0=time.time()
#     for ds in datasets:
#         name, path = ds['name'], ds['path']
#         files = sorted(glob.glob(f"{path}/*.root"))
#         if args.maxFiles>0: files = files[:args.maxFiles]
#         if not files:
#             log.warning("%s   no files – skipped", name); continue
#         log.info("Dataset %s: %d files", name, len(files))

#         df = ROOT.RDataFrame("events", files)
#         _,H,nEvPtr,occPtr = my_graph(df, name)
#         graphs.extend(H + [nEvPtr, occPtr])
#         per_ds.append((name,H,nEvPtr,occPtr))

#     if not graphs:
#         log.error("Nothing to process – exit"); return

#     log.info("All graphs built – starting event loop …")
#     ROOT.RDF.RunGraphs(graphs)
#     log.info("Event loop finished in %.1f s", time.time()-t0)

#     out_dir = pathlib.Path(args.outputDir); out_dir.mkdir(exist_ok=True)

#     for name,H,nEvPtr,occPtr in per_ds:
#         nEv   = int(nEvPtr.GetValue())
#         occ_h = occPtr.GetValue()
#         raw_hits = occ_h.Integral()

#         # save ROOT file
#         fout = ROOT.TFile(str(out_dir/f"{name}.root"),"RECREATE")
#         for hptr in H:
#             h = hptr.GetValue(); h.Scale(1.0/nEv); h.Write()
#         occ_h.Write()
#         # ints / floats stored with TParameter
#         ROOT.TParameter(int)   ("nEvents",    nEv     ).Write()
#         ROOT.TParameter(int)   ("occ_hits",   int(raw_hits)).Write()

#         vals = [occ_h.GetBinContent(i)/nEv
#                 for i in range(1,occ_h.GetNbinsX()+1) if occ_h.GetBinContent(i)]
#         sf,cl = 3,5
#         occ_max = max(vals)*sf*cl if vals else 0.
#         occ_avg = (sum(vals)/len(vals))*sf*cl if vals else 0.
#         ROOT.TParameter(float) ("occ_max",    occ_max).Write()
#         ROOT.TParameter(float) ("occ_avg",    occ_avg).Write()
#         fout.Close()

#         log.info("[%s] occupancy hist integral %d (%.2f per event)",
#                  name, raw_hits, raw_hits/nEv)

# # plotting & stats
# def make_plots(rfile, tag, *, logz=False):
#     f = ROOT.TFile(rfile)
#     base = pathlib.Path(rfile).with_suffix('')
#     dH = base.parent / f"hits_{tag}";   dH.mkdir(exist_ok=True)
#     dE = base.parent / f"energy_{tag}"; dE.mkdir(exist_ok=True)
#     dM = base.parent / f"mc_{tag}";     dM.mkdir(exist_ok=True)

#     for n,lbl in [("z","z (mm) / event"),("theta","theta (deg) / event"),
#                   ("phi","phi (deg) / event"),("z_layer0","z L0 / event"),
#                   ("theta_layer0","theta L0 / event"),
#                   ("phi_layer0","phi L0 / event")]:
#         plot_hist(f.Get(n), f"{dH}/{n}.png", lbl)
#     for n,lbl in [("z_phi","z vs phi / event"),
#                   ("theta_phi","theta vs phi / event"),
#                   ("z_phi_layer0","z vs phi L0 / event"),
#                   ("theta_phi_layer0","theta vs phi L0 / event")]:
#         plot_2d(f.Get(n), f"{dH}/{n}.png", lbl, logz)
#     for i,_ in enumerate(layer_radii,1):
#         ep=f.Get("energy").ProjectionY(f"E{i}",i,i)
#         dx=f.Get("dedx").  ProjectionY(f"dE{i}",i,i)
#         plot_hist(ep,f"{dE}/energy_L{i}.png",f"Energy L{i} / event")
#         plot_hist(dx,f"{dE}/dedx_L{i}.png",  f"dE/dx L{i} / event")
#     plot_2d(f.Get("mc_z_r"),      f"{dM}/mc_z_r.png",
#             "MC z vs r / event",logz)
#     plot_2d(f.Get("mc_theta_phi"),f"{dM}/mc_theta_phi.png",
#             "MC theta vs phi / event",logz)
#     plot_hist(f.Get("mc_energy"), f"{dM}/mc_energy.png",
#               "MC energy / event")

# def occupancy_stats(rfile, tag):
#     f      = ROOT.TFile(rfile)
#     occMax = f.Get("occ_max")
#     occAvg = f.Get("occ_avg")
#     if occMax and occAvg:          # preferred (fast) path
#         print(f"[{tag}]  Max occ = {occMax.GetVal():.3f}  "
#               f"Avg occ = {occAvg.GetVal():.3f}")
#         return
#     # fallback: compute on the fly
#     h   = f.Get("cell_counts")
#     nEv = f.Get("nEvents").GetVal()
#     vals=[h.GetBinContent(i)/nEv for i in range(1,h.GetNbinsX()+1)
#           if h.GetBinContent(i)]
#     sf,cl=3,5
#     print(f"[{tag}]  Max occ = {max(vals)*sf*cl:.3f}  "
#           f"Avg occ = {sum(vals)/len(vals)*sf*cl:.3f}")


# DATASETS = [
#     {"name":"256",
#      "path":"/ceph/submit/data/user/k/kudela/beam_backgrounds/CLD_o2_v05/FCCee_Z_4IP_04may23_FCCee_Z256"},
#     {"name":"128",
#      "path":"/ceph/submit/data/user/k/kudela/beam_backgrounds/CLD_o2_v05/FCCee_Z_4IP_04may23_FCCee_Z128"},
# ]


# def main():
#     out_dir = pathlib.Path(args.outputDir); out_dir.mkdir(exist_ok=True)

#     if args.calculate:
#         build_and_run(DATASETS)

#     for ds in DATASETS:
#         rootfile = out_dir/f"{ds['name']}.root"
#         if args.plots:
#             make_plots(str(rootfile), ds['name'], logz=args.log)
#         if args.stats:
#             occupancy_stats(str(rootfile), ds['name'])

# if __name__ == "__main__":
#     main()
