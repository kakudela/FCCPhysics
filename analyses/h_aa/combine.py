import sys,array,ROOT,math,os,copy
import argparse
import numpy as np

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

parser = argparse.ArgumentParser()
parser.add_argument("--run", help="Run combine", action='store_true')
parser.add_argument("--ecm", type=int, help="Center-of-mass energy", default=240)
args = parser.parse_args()

def getHist(hName, procs, rebin=1):
    hist = None
    for proc in procs:
        fIn = ROOT.TFile(f"{inputDir}/{proc}.root")
        h = fIn.Get(hName)
        h.SetDirectory(0)
        if hist == None:
            hist = h
        else:
            hist.Add(h)
        fIn.Close()
    hist.Rebin(rebin)
    return hist


def smooth_histogram_gaussian_kernel(hist, sigma, fOut):

    n_bins = hist.GetNbinsX()
    smoothed_values = []

    # Loop over all bins
    for i in range(1, n_bins + 1):
        bin_center = hist.GetBinCenter(i)
        kernel_sum = 0.0
        value_sum = 0.0

        # Apply Gaussian kernel to surrounding bins
        for j in range(1, n_bins + 1):
            neighbor_center = hist.GetBinCenter(j)
            weight = np.exp(-0.5 * ((bin_center - neighbor_center) / sigma) ** 2)
            value_sum += hist.GetBinContent(j) * weight
            kernel_sum += weight

        # Normalize the value by the kernel sum
        smoothed_values.append(value_sum / kernel_sum)

    smoothed_hist = hist.Clone(f"{hist.GetName()}_smoothed")
    smoothed_hist.Reset()  # Clear bin contents
    for i, value in enumerate(smoothed_values, start=1):
        smoothed_hist.SetBinContent(i, value)


    canvas = ROOT.TCanvas("", "", 800, 800)

    # Draw the original histogram
    hist.SetLineColor(ROOT.kBlue)
    hist.SetLineWidth(2)
    hist.Draw("HIST")

    smoothed_hist.SetLineColor(ROOT.kRed)
    smoothed_hist.SetLineWidth(2)
    smoothed_hist.Draw("HIST SAME")

    legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
    legend.AddEntry(hist, "Original", "l")
    legend.AddEntry(smoothed_hist, "Smoothed", "l")
    legend.Draw()

    canvas.SaveAs(f"{plot_dir}/{fOut}.png")
    canvas.SaveAs(f"{plot_dir}/{fOut}.pdf")


    return smoothed_hist

if __name__ == "__main__":

    ecm = args.ecm

    inputDir = f"output/h_aa/histmaker/ecm{ecm}/smearIDEA_S10pct"
    outDir = f"output/h_aa/combine/ecm{ecm}/"
    plot_dir = f"/home/submit/jaeyserm/public_html/fccee/h_aa/combine_smoothing_ecm{ecm}/"


    if ecm == 240:
        sigs = ['wzp6_ee_nunuH_Haa_ecm240', 'wzp6_ee_eeH_Haa_ecm240', 'wzp6_ee_tautauH_Haa_ecm240', 'wzp6_ee_ccH_Haa_ecm240', 'wzp6_ee_bbH_Haa_ecm240', 'wzp6_ee_qqH_Haa_ecm240', 'wzp6_ee_ssH_Haa_ecm240', 'wzp6_ee_mumuH_Haa_ecm240']
        bkgs = ['wz3p6_ee_uu_ecm240', 'wz3p6_ee_dd_ecm240', 'wz3p6_ee_cc_ecm240', 'wz3p6_ee_ss_ecm240', 'wz3p6_ee_bb_ecm240', 'wz3p6_ee_tautau_ecm240', 'wz3p6_ee_mumu_ecm240', 'wz3p6_ee_nunu_ecm240', 'wz3p6_ee_gammagamma_ecm240', 'wz3p6_ee_ee_Mee_30_150_ecm240']

    if ecm == 365:
        sigs = ['wzp6_ee_eeH_Haa_ecm365', 'wzp6_ee_tautauH_Haa_ecm365', 'wzp6_ee_ccH_Haa_ecm365', 'wzp6_ee_bbH_Haa_ecm365', 'wzp6_ee_qqH_Haa_ecm365', 'wzp6_ee_ssH_Haa_ecm365', 'wzp6_ee_mumuH_Haa_ecm365'] # remove nunuH as treated separately
        bkgs = ['wz3p6_ee_uu_ecm365', 'wz3p6_ee_dd_ecm365', 'wz3p6_ee_cc_ecm365', 'wz3p6_ee_ss_ecm365', 'wz3p6_ee_bb_ecm365', 'wz3p6_ee_tautau_ecm365', 'wz3p6_ee_mumu_ecm365', 'wz3p6_ee_nunu_ecm365', 'wz3p6_ee_gammagamma_ecm365', 'wz3p6_ee_ee_Mee_30_150_ecm365']

    cats = ['zmumu', 'zee', 'znunu', 'zqq']
    if ecm == 365:
        cats.append('vbf')
    for cat in cats:
        print("**************************", cat)
        if ecm == 240:
            rebin = 2
            sigma = 1.25
            hName = f"{cat}_haa_m"
            h_sig = getHist(hName, sigs, rebin=rebin)
            h_bkg = getHist(hName, bkgs, rebin=rebin)
            #if cat == "ee":
            #    h_bkg = getHist(hName.replace("ee", "mumu"), bkgs, rebin=rebin) # electron as muon background
            #    h_bkg.Scale(5.3) # ratio Z->ee / Z->mumu at 365 (8.31E+03/1.57E+03)

            h_bkg = smooth_histogram_gaussian_kernel(h_bkg, sigma, f"{cat}_bkg")

            h_sig.SetName(f"{cat}_sig")
            h_bkg.SetName(f"{cat}_bkg")

            fOut = ROOT.TFile(f"{outDir}/datacard_{cat}.root", "RECREATE")
            h_sig.Write()
            h_bkg.Write()
            h_data = h_sig.Clone(f"{cat}_data")
            h_data.Add(h_bkg)
            h_data.Write()
            fOut.Close()

            dc = ""
            dc += "imax *\n"
            dc += "jmax *\n"
            dc += "kmax *\n"
            dc += "####################\n"
            dc += f"shapes *        * datacard_{cat}.root $CHANNEL_$PROCESS\n"
            dc += f"shapes data_obs * datacard_{cat}.root $CHANNEL_data\n"
            dc += "####################\n"
            dc += f"bin          {cat}\n"
            dc += "observation  -1\n"
            dc += "####################\n"
            dc += f"bin          {cat}       {cat}\n"
            dc += "process      sig         bkg\n"
            dc += "process      0           1\n"
            dc += "rate         -1          -1\n"
            dc += "####################\n"
            #dc += "dummy lnN    -          1.00001\n"
            dc += "bkg lnN      -           1.05\n"

        if ecm == 365:
            rebin = 2
            sigma = 1.75
            hName = f"{cat}_haa_m"
            h_sig_zh = getHist(hName, sigs, rebin=rebin)
            h_sig_numunumuH = getHist(hName, ['wzp6_ee_numunumuH_Haa_ecm365'], rebin=rebin)
            h_sig_numunumuH.Scale(3)
            h_sig_zh.Add(h_sig_numunumuH)

            h_sig_vbf = getHist(hName, ['wzp6_ee_nuenueH_Haa_ecm365'], rebin=rebin)
            h_sig_numunumuH = getHist(hName, ['wzp6_ee_numunumuH_Haa_ecm365'], rebin=rebin)
            print(h_sig_vbf.Integral())
            print(h_sig_numunumuH.Integral())
            h_sig_vbf.Add(h_sig_numunumuH, -1)
            print(h_sig_vbf.Integral())

            for k in range(1, h_sig_vbf.GetNbinsX()+1):
                if h_sig_vbf.GetBinContent(k) <= 0:
                    h_sig_vbf.SetBinContent(k, 1e-5)

            h_bkg = getHist(hName, bkgs, rebin=rebin)
            h_bkg = smooth_histogram_gaussian_kernel(h_bkg, sigma, f"{cat}_bkg")
            #h_bkg.Scale(0.0001)
            #h_sig_zh.Scale(0.0001)

            h_sig_zh.SetName(f"{cat}_sig_zh")
            h_sig_vbf.SetName(f"{cat}_sig_vbf")
            h_bkg.SetName(f"{cat}_bkg")

            fOut = ROOT.TFile(f"{outDir}/datacard_{cat}.root", "RECREATE")
            h_sig_zh.Write()
            h_sig_vbf.Write()
            h_bkg.Write()
            h_data = h_sig_zh.Clone(f"{cat}_data")
            h_data.Add(h_sig_vbf)
            h_data.Add(h_bkg)
            h_data.Write()
            fOut.Close()


            dc = ""
            dc += "imax *\n"
            dc += "jmax *\n"
            dc += "kmax *\n"
            dc += "####################\n"
            dc += f"shapes *        * datacard_{cat}.root $CHANNEL_$PROCESS\n"
            dc += f"shapes data_obs * datacard_{cat}.root $CHANNEL_data\n"
            dc += "####################\n"
            dc += f"bin          {cat}\n"
            dc += "observation  -1\n"
            dc += "####################\n"
            dc += f"bin          {cat}      {cat}      {cat}\n"
            dc += "process      sig_zh      sig_vbf     bkg\n"
            dc += "process      0           -2          1\n"
            dc += "rate         -1          -1          -1\n"
            dc += "####################\n"
            #dc += "dummy lnN    -     -     1.00001\n"
            dc += "bkg lnN      -           -           1.05\n"

        f = open(f"{outDir}/datacard_{cat}.txt", 'w')
        f.write(dc)
        f.close()

        print(dc)

        if args.run:
            cmd = f"singularity exec --bind /work:/work /work/submit/jaeyserm/software/docker/cmssw_cc7.sif bash -c 'PYTHONPATH=''; dir=$(pwd); cd /work/submit/jaeyserm/wmass/CMSSW_10_6_19_patch2/src; source /cvmfs/cms.cern.ch/cmsset_default.sh; cmsenv; cd $dir/{outDir}; text2hdf5.py --X-allow-no-background datacard_{cat}.txt -o ws_{cat}.hdf5; combinetf.py ws_{cat}.hdf5 -o fit_output_{cat}.root -t 0  --expectSignal=1'"
            os.system(cmd)
    #quit()

    if args.run:
        cards = ' '.join([f'datacard_{x}.txt' for x in cats])
        cmd = f"singularity exec --bind /work:/work /work/submit/jaeyserm/software/docker/cmssw_cc7.sif bash -c 'PYTHONPATH=''; dir=$(pwd); cd /work/submit/jaeyserm/wmass/CMSSW_10_6_19_patch2/src; source /cvmfs/cms.cern.ch/cmsset_default.sh; cmsenv; cd $dir/{outDir}; combineCards.py {cards} > datacard_combined.txt'"
        os.system(cmd)

        cmd = f"singularity exec --bind /work:/work /work/submit/jaeyserm/software/docker/cmssw_cc7.sif bash -c 'PYTHONPATH=''; dir=$(pwd); cd /work/submit/jaeyserm/wmass/CMSSW_10_6_19_patch2/src; source /cvmfs/cms.cern.ch/cmsset_default.sh; cmsenv; cd $dir/{outDir}; text2hdf5.py --X-allow-no-background datacard_combined.txt -o ws_combined.hdf5; combinetf.py ws_combined.hdf5 -o fit_output_combined.root -t 0  --expectSignal=1 '" # --binByBinStat
        os.system(cmd)

        for cat in cats+['combined']:
            fIn = ROOT.TFile(f"{outDir}/fit_output_{cat}.root")
            tree = fIn.Get("fitresults")
            tree.GetEntry(0)
            status = tree.status
            errstatus = tree.errstatus
            
            #mu = tree.sig_vbf_mu
            #mu_err = tree.sig_vbf_mu_err*100.
            #print(f"{cat}\tzh={mu_err:.3f} {status}")
            #continue
            if ecm == 365:
                mu_zh = tree.sig_zh_mu
                mu_zh_err = tree.sig_zh_mu_err*100.
                mu_vbf = tree.sig_vbf_mu
                mu_vbf_err = tree.sig_vbf_mu_err*100.
                print(f"{cat}\tzh={mu_zh_err:.3f}\tvbf={mu_vbf_err:.3f} {status}")
            else:
                mu_zh = tree.sig_mu
                mu_zh_err = tree.sig_mu_err*100.
                print(f"{cat}\tzh={mu_zh_err:.3f} {status}")

