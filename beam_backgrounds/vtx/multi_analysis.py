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
        # MC 
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
