import ROOT, os, argparse, logging
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

"""
This script processes MC data from ROOT files, generating heat-maps and 1D energy distributions.
The outdirectory is mc_analysis.
"""

ap = argparse.ArgumentParser()
ap.add_argument('--tags',   nargs='+', default=['256','128'], help='dataset tags (expects output/<tag>.root)')
ap.add_argument('--outDir', default='mc_analysis')
ap.add_argument('--log',    action='store_true', help='draw heat-maps with log colour scale')
args = ap.parse_args()

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
log = logging.getLogger('mc')

os.makedirs(f'{args.outDir}/single', exist_ok=True)

# helpers
def canvas():
    c = ROOT.TCanvas('', '', 800, 650)
    c.SetLeftMargin(0.12); c.SetRightMargin(0.15)
    return c

def draw_heatmap(h2, fname, logz):
    c = canvas()
    if logz: c.SetLogz()
    h2.Draw('COLZ')
    c.SaveAs(fname); c.Close()

def draw_1d(h1, xlab, title, fname):
    c = ROOT.TCanvas('', '', 800, 600)
    h1.SetLineWidth(2); h1.SetLineColor(ROOT.kBlue+1)
    h1.SetTitle(title); h1.GetXaxis().SetTitle(xlab)
    h1.GetYaxis().SetTitle('MC particles / event')
    h1.Draw('HIST')
    c.SaveAs(fname); c.Close()

# loop
for tag in args.tags:
    fpath = f'output/{tag}.root'
    f = ROOT.TFile(fpath)
    if not f or f.IsZombie():
        log.warning(f'skip {tag}: cannot open {fpath}')
        continue

    nEv = f.Get('nEvents').GetVal()
    if nEv == 0:
        log.warning(f'skip {tag}: nEvents == 0 in {fpath}')
        f.Close(); continue
    scale = 1./nEv

    # grab, clone, detach so file can be closed
    hist = {}
    for k in ('mc_z_r', 'mc_z_phi', 'mc_theta_phi', 'mc_E'):
        obj = f.Get(k)
        if not obj: log.warning(f'{k} not found in {tag}'); continue
        h = obj.Clone(f'{k}_{tag}'); h.SetDirectory(0); h.Scale(scale)
        hist[k] = h
    f.Close()

    # heat-maps
    if 'mc_z_r' in hist:
        h = hist['mc_z_r']
        h.SetTitle(f'z vs r ({tag})')
        h.GetXaxis().SetRangeUser(-60, 60)
        h.GetYaxis().SetRangeUser(0, 50)
        h.GetXaxis().SetTitle('z (mm)')
        h.GetYaxis().SetTitle('r (mm)')
        draw_heatmap(h, f'{args.outDir}/single/z_r_{tag}.png', args.log)

    if 'mc_z_phi' in hist:
        h = hist['mc_z_phi']
        h.SetTitle(f'z vs phi ({tag})')
        h.GetXaxis().SetRangeUser(-60, 60)
        h.GetXaxis().SetTitle('z (mm)')
        h.GetYaxis().SetTitle('phi (deg)')
        draw_heatmap(h, f'{args.outDir}/single/z_phi_{tag}.png', args.log)

    if 'mc_theta_phi' in hist:
        h = hist['mc_theta_phi']
        h.SetTitle(f'theta vs phi ({tag})')
        h.GetXaxis().SetTitle('theta (deg)')
        h.GetYaxis().SetTitle('phi (deg)')
        draw_heatmap(h, f'{args.outDir}/single/theta_phi_{tag}.png', args.log)

    # energy 1-D
    if 'mc_E' in hist:
        hist['mc_E'].GetXaxis().SetRangeUser(0, 15)
        draw_1d(hist['mc_E'], 'energy (GeV)',
                f'Energy distribution ({tag})',
                f'{args.outDir}/single/energy_{tag}.png')

log.info('MC plotting finished.')
