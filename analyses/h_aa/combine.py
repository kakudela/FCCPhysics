import sys,array,ROOT,math,os,copy
import argparse
import numpy as np

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

parser = argparse.ArgumentParser()
parser.add_argument("--cat", type=str, help="Category (mumu, ee, nunu, qq)", default="mumu")
parser.add_argument("--run", help="Run combine", action='store_true')
parser.add_argument("--plot", help="Plot", action='store_true')
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

    inputDir = "output/h_aa/histmaker/"
    outDir = "output/h_aa/combine/"
    plot_dir = "/home/submit/jaeyserm/public_html/fccee/h_aa/combine_smoothing/"
    rebin = 2
    sigma = 1

    sigs = ['wzp6_ee_nunuH_Haa_ecm240', 'wzp6_ee_eeH_Haa_ecm240', 'wzp6_ee_tautauH_Haa_ecm240', 'wzp6_ee_ccH_Haa_ecm240', 'wzp6_ee_bbH_Haa_ecm240', 'wzp6_ee_qqH_Haa_ecm240', 'wzp6_ee_ssH_Haa_ecm240', 'wzp6_ee_mumuH_Haa_ecm240']
    bkgs = ['kkmcee_ee_uu_ecm240', 'kkmcee_ee_dd_ecm240', 'kkmcee_ee_cc_ecm240', 'kkmcee_ee_ss_ecm240', 'kkmcee_ee_bb_ecm240', 'kkmcee_ee_tautau_ecm240', 'kkmcee_ee_mumu_ecm240', 'kkmcee_ee_nuenue_ecm240', 'kkmcee_ee_numunumu_ecm240', 'kkmcee_ee_nutaunutau_ecm240', 'wzp6_ee_gammagamma_ecm240']
    bkgs = ['wz3p6_ee_uu_ecm240', 'wz3p6_ee_dd_ecm240', 'wz3p6_ee_cc_ecm240', 'wz3p6_ee_ss_ecm240', 'wz3p6_ee_bb_ecm240', 'wz3p6_ee_tautau_ecm240', 'wz3p6_ee_mumu_ecm240', 'wz3p6_ee_nunu_ecm240', 'wzp6_ee_gammagamma_ecm240', 'wz3p6_ee_ee_Mee_30p_ecm240']


    for cat in ['mumu', 'ee', 'nunu', 'qq']:
        hName = f"z{cat}_haa_m"
        h_sig = getHist(hName, sigs, rebin=rebin)
        h_bkg = getHist(hName, bkgs, rebin=rebin)
        if cat == "ee":
            h_bkg = getHist(hName.replace("ee", "mumu"), bkgs, rebin=rebin) # electron as muon background
            h_bkg.Scale(5.3) # ratio Z->ee / Z->mumu at 240 (8.31E+03/1.57E+03)

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
        dc += "dummy lnN    1.00001     1.00001\n"
        dc += "bkg lnN      -           1.50\n"
        #dc += "* autoMCStats 0 0 1\n"


        f = open(f"{outDir}/datacard_{cat}.txt", 'w')
        f.write(dc)
        f.close()

        print(dc)

        if args.run:
            #cmd = f"singularity exec --bind /work:/work /work/submit/jaeyserm/software/docker/combine-standalone_v9.2.1.sif bash -c 'cd {outDir}; text2workspace.py datacard_{cat}.txt -o ws_{cat}.root; combine -M MultiDimFit -v 10 --rMin 0.5 --rMax 1.5 --setParameters r=1 ws_{cat}.root'"
            cmd = f"singularity exec --bind /work:/work /work/submit/jaeyserm/software/docker/combine-standalone_v9.2.1.sif bash -c 'cd {outDir}; text2workspace.py datacard_{cat}.txt -o ws_{cat}.root; combine -M FitDiagnostics -t -1 --setParameters r=1 ws_{cat}.root -n {cat} --cminDefaultMinimizerStrategy 0'"
            os.system(cmd)

    if args.run:
        cmd = f"singularity exec --bind /work:/work /work/submit/jaeyserm/software/docker/combine-standalone_v9.2.1.sif bash -c 'cd {outDir}; combineCards.py datacard_mumu.txt datacard_ee.txt datacard_nunu.txt datacard_qq.txt > datacard_combined.txt'"
        os.system(cmd)

        cmd = f"singularity exec --bind /work:/work /work/submit/jaeyserm/software/docker/combine-standalone_v9.2.1.sif bash -c 'cd {outDir}; text2workspace.py datacard_combined.txt -o ws_combined.root; combine -M FitDiagnostics -t -1 --setParameters r=1 ws_combined.root -n combined --cminDefaultMinimizerStrategy 0'"
        os.system(cmd)

        ## extract results
        for cat in ['mumu', 'ee', 'nunu', 'qq', 'combined']:
            fIn = ROOT.TFile(f"{outDir}/fitDiagnostics{cat}.root")
            t = fIn.Get("tree_fit_sb")
            t.GetEntry(0)
            err = t.rErr
            status = t.fit_status
            print(f"{cat}\t{err:.3f} {status}")