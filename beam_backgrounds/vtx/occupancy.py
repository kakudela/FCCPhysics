import ROOT, argparse, os, glob, logging
ROOT.gROOT.SetBatch(True)

"""
This script processes occupancy data from ROOT files (multi_analysis.py), calculating maximum and average occupancy
"""

ap = argparse.ArgumentParser()
ap.add_argument('--inputDir', default='output', help='directory that holds <tag>.root files from multi_analysis')
ap.add_argument('--tags', nargs='+', default=None, help='dataset tags; if omitted, use every *.root in inputDir')
ap.add_argument('--safety',  type=float, default=3, help='safety factor')
ap.add_argument('--cluster', type=float, default=5, help='cluster size (pixels per module)')
args = ap.parse_args()

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
log = logging.getLogger('occ')

# find the root files we will inspect
if args.tags:
    files = [f"{args.inputDir}/{tag}.root" for tag in args.tags]
else:
    files = glob.glob(f"{args.inputDir}/*.root")

if not files:
    log.error('No ROOT files found – nothing to do.')
    exit()

for path in files:
    tag = os.path.splitext(os.path.basename(path))[0]
    f = ROOT.TFile(path)
    if not f or f.IsZombie():
        log.warning(f"skip {tag}: cannot open {path}")
        continue

    h = f.Get('cell_counts')
    par = f.Get('nEvents')
    if not h or not par:
        log.warning(f"skip {tag}: required objects not found")
        f.Close(); continue

    nEv = par.GetVal()
    vals = [h.GetBinContent(i)/nEv for i in range(1, h.GetNbinsX()+1)
            if h.GetBinContent(i)]
    if not vals:
        log.warning(f"skip {tag}: empty occupancy histogram")
        f.Close(); continue

    max_occ = max(vals) * args.safety * args.cluster
    avg_occ = (sum(vals)/len(vals)) * args.safety * args.cluster

    print(f"{tag:>8}  Max occupancy: {max_occ:.3f}   Avg occupancy: {avg_occ:.3f}")
    f.Close()
