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
    #hist.Rebin(rebin)
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


def unroll_2d(hist):

    # Get the number of bins in X and Y
    n_bins_x = hist.GetNbinsX()
    n_bins_y = hist.GetNbinsY()

    # Create a 1D histogram to hold the unrolled data
    n_bins_1d = n_bins_x * n_bins_y
    h1 = ROOT.TH1D("h1", "1D Unrolled Histogram", n_bins_1d, 0.5, n_bins_1d + 0.5)

    # Loop over all bins in the 2D histogram
    for bin_x in range(1, n_bins_x + 1):  # Bin indexing starts at 1
        for bin_y in range(1, n_bins_y + 1):
            # Calculate the global bin number for the 1D histogram
            bin_1d = (bin_y - 1) * n_bins_x + bin_x

            # Get the content and error from the 2D histogram
            content = hist.GetBinContent(bin_x, bin_y)
            error = hist.GetBinError(bin_x, bin_y)

            # Set the content and error in the 1D histogram
            h1.SetBinContent(bin_1d, content)
            h1.SetBinError(bin_1d, error)
    return h1

if __name__ == "__main__":

    inputDir = "output/h_za/histmaker/"
    outDir = "output/h_za/combine/"
    plot_dir = "/home/submit/jaeyserm/public_html/fccee/h_za/combine_smoothing/"
    rebin = 1
    sigma = 1

    sigs = ['wzp6_ee_nunuH_HZa_ecm240', 'wzp6_ee_eeH_HZa_ecm240', 'wzp6_ee_tautauH_HZa_ecm240', 'wzp6_ee_ccH_HZa_ecm240', 'wzp6_ee_bbH_HZa_ecm240', 'wzp6_ee_qqH_HZa_ecm240', 'wzp6_ee_ssH_HZa_ecm240', 'wzp6_ee_mumuH_HZa_ecm240']
    bkgs = ['kkmcee_ee_uu_ecm240', 'kkmcee_ee_dd_ecm240', 'kkmcee_ee_cc_ecm240', 'kkmcee_ee_ss_ecm240', 'kkmcee_ee_bb_ecm240', 'kkmcee_ee_tautau_ecm240', 'kkmcee_ee_mumu_ecm240', 'kkmcee_ee_nuenue_ecm240', 'kkmcee_ee_numunumu_ecm240', 'kkmcee_ee_nutaunutau_ecm240', 'wzp6_ee_gammagamma_ecm240']
    bkgs = ['p8_ee_Zqq_ecm240']
    #bkgs = ['kkmcee_ee_uu_ecm240', 'kkmcee_ee_dd_ecm240', 'kkmcee_ee_cc_ecm240', 'kkmcee_ee_ss_ecm240', 'kkmcee_ee_bb_ecm240']
    bkgs = ['wz3p6_ee_uu_ecm240', 'wz3p6_ee_dd_ecm240', 'wz3p6_ee_cc_ecm240', 'wz3p6_ee_ss_ecm240', 'wz3p6_ee_bb_ecm240', 'wz3p6_ee_tautau_ecm240', 'wz3p6_ee_mumu_ecm240', 'wz3p6_ee_nunu_ecm240', 'wzp6_ee_gammagamma_ecm240']

    hists = ['final_histo'] # 
    hists = ["vvqqqqa_final", "vvqqvva_final"]

    for hName in hists:
        h_sig = getHist(hName, sigs, rebin=rebin)
        h_bkg = getHist(hName, bkgs, rebin=rebin)

        #h_sig = unroll_2d(h_sig)
        #h_bkg = unroll_2d(h_bkg)
        #h_sig = smooth_histogram_gaussian_kernel(h_sig, sigma, f"{hName}_sig")
        h_bkg = smooth_histogram_gaussian_kernel(h_bkg, sigma, f"{hName}_bkg")

        h_sig.SetName(f"{hName}_sig")
        h_bkg.SetName(f"{hName}_bkg")

        fOut = ROOT.TFile(f"{outDir}/datacard_{hName}.root", "RECREATE")
        h_sig.Write()
        h_bkg.Write()
        h_data = h_sig.Clone(f"{hName}_data")
        h_data.Add(h_bkg)
        h_data.Write()
        fOut.Close()

        dc = ""
        dc += "imax *\n"
        dc += "jmax *\n"
        dc += "kmax *\n"
        dc += "####################\n"
        dc += f"shapes *        * datacard_{hName}.root $CHANNEL_$PROCESS\n"
        dc += f"shapes data_obs * datacard_{hName}.root $CHANNEL_data\n"
        dc += "####################\n"
        dc += f"bin          {hName}\n"
        dc += "observation  -1\n"
        dc += "####################\n"
        dc += f"bin          {hName}       {hName}\n"
        dc += "process      sig         bkg\n"
        dc += "process      0           1\n"
        dc += "rate         -1          -1\n"
        dc += "####################\n"
        dc += "dummy lnN    1.00001     1.00001\n"
        #dc += "bkg lnN      -           1.50\n"
        #dc += "* autoMCStats 0 0 1\n"


        f = open(f"{outDir}/datacard_{hName}.txt", 'w')
        f.write(dc)
        f.close()

        print(dc)

        if args.run:
            #cmd = f"singularity exec --bind /work:/work /work/submit/jaeyserm/software/docker/combine-standalone_v9.2.1.sif bash -c 'cd {outDir}; text2workspace.py datacard_{cat}.txt -o ws_{cat}.root; combine -M MultiDimFit -v 10 --rMin 0.5 --rMax 1.5 --setParameters r=1 ws_{cat}.root'"
            cmd = f"singularity exec --bind /work:/work /work/submit/jaeyserm/software/docker/combine-standalone_v9.2.1.sif bash -c 'cd {outDir}; text2workspace.py datacard_{hName}.txt -o ws_{hName}.root; combine -M FitDiagnostics -t -1 --setParameters r=1 ws_{hName}.root -n {hName} --cminDefaultMinimizerStrategy 0'"
            os.system(cmd)

    if args.run:
        cards = ' '.join([f'datacard_{x}.txt' for x in hists])
        cmd = f"singularity exec --bind /work:/work /work/submit/jaeyserm/software/docker/combine-standalone_v9.2.1.sif bash -c 'cd {outDir}; combineCards.py {cards} > datacard_combined.txt'"
        os.system(cmd)

        cmd = f"singularity exec --bind /work:/work /work/submit/jaeyserm/software/docker/combine-standalone_v9.2.1.sif bash -c 'cd {outDir}; text2workspace.py datacard_combined.txt -o ws_combined.root; combine -M FitDiagnostics -t -1 --setParameters r=1 ws_combined.root -n combined --cminDefaultMinimizerStrategy 0'"
        os.system(cmd)

        ## extract results
        for hist in hists + ['combined']:
            fIn = ROOT.TFile(f"{outDir}/fitDiagnostics{hist}.root")
            t = fIn.Get("tree_fit_sb")
            t.GetEntry(0)
            err = t.rErr
            status = t.fit_status
            print(f"{hist}\t{err:.3f} {status}")

