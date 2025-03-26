import sys,array,ROOT,math,os,copy
import argparse
import numpy as np

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

parser = argparse.ArgumentParser()
parser.add_argument("--ecm", type=int, help="Center-of-mass energy", default=240)
parser.add_argument("--run", help="Run combine", action='store_true')
args = parser.parse_args()

def getHistQQ(hName, procs, rebin=1):
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

def removeNegativeBins(hist):
    totNeg, tot = 0., 0.
    if "TH1" in hist.ClassName():
        nbinsX = hist.GetNbinsX()
        for x in range(1, nbinsX + 1):
            content = hist.GetBinContent(x)
            tot += content
            if content < 0:
                totNeg += content
                hist.SetBinContent(x, 0)
                hist.SetBinError(x, 0)
    elif "TH2" in hist.ClassName():
        pass

    elif "TH3" in hist.ClassName():
        nbinsX = hist.GetNbinsX()
        nbinsY = hist.GetNbinsY()
        nbinsZ = hist.GetNbinsZ()
        for x in range(1, nbinsX + 1):
            for y in range(1, nbinsY + 1):
                for z in range(1, nbinsZ + 1):
                    content = hist.GetBinContent(x, y, z)
                    error = hist.GetBinError(x, y, z)  # Retrieve bin error
                    tot += content
                    if content < 0:
                        #print("WARNING: NEGATIVE BIN CONTENT", content, hist.GetName())
                        totNeg += content
                        hist.SetBinContent(x, y, z, 0)
                        hist.SetBinError(x, y, z, 0)
    if totNeg != 0:
        print(f"WARNING: TOTAL {tot}, NEGATIVE {totNeg}, FRACTION {totNeg/tot}")
    return hist


def getHist(hName, procs, rebin=1, isVBF=False):
    hist = None
    for proc in procs:
        #isVBF = False
        if 'nuenueVBF' in proc:
            proc = proc.replace('VBF', '')
            isVBF = True
        fIn = ROOT.TFile(f"{inputDir}/{proc}.root")
        h = fIn.Get(hName)
        h.SetDirectory(0)
        fIn.Close()
        if 'nuenue' in proc and isVBF:
            h_numunumu = getHist(hName, [proc.replace('nuenue', 'numunumu')], isVBF=True)
            h.Add(h_numunumu, -1)
        if "numunumuH" in proc and not isVBF:
            h.Scale(3.)
        if hist == None:
            hist = h
        else:
            hist.Add(h)
        
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

def unroll(hist, rebin=1):
    if "TH1" in hist.ClassName():
        hist = hist.Rebin(rebin)
        hist = removeNegativeBins(hist)
        return hist

    elif "TH2" in hist.ClassName():
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
                #print(hist.GetXaxis().GetBinCenter(bin_x), hist.GetYaxis().GetBinCenter(bin_y), content, error)
        h1 = removeNegativeBins(h1)
        return h1


    elif "TH3" in hist.ClassName():
        # Get binning information
        nbinsX = hist.GetNbinsX()
        nbinsY = hist.GetNbinsY()
        nbinsZ = hist.GetNbinsZ()
        nbins1D = nbinsX * nbinsY * nbinsZ

        # Create a 1D histogram with the correct number of bins
        h1 = ROOT.TH1D("h1_unrolled", "Unrolled 3D Histogram", nbins1D, 0, nbins1D)

        # Fill the 1D histogram by unrolling the 3D histogram
        bin1D = 1  # ROOT bins are 1-based
        for x in range(1, nbinsX + 1):
            for y in range(1, nbinsY + 1):
                for z in range(1, nbinsZ + 1):
                    content = hist.GetBinContent(x, y, z)
                    error = hist.GetBinError(x, y, z)  # Retrieve bin error
                    
                    if content < 0:
                        print("WARNING NEGATIVE CONTENT", content, hist.GetName())
                        content = 0
                        error = 0

                    h1.SetBinContent(bin1D, content)
                    h1.SetBinError(bin1D, error)  # Set the bin error
                    bin1D += 1  # Increment the 1D bin index
                    
        return h1

    else:
        return hist


if __name__ == "__main__":

    ecm = args.ecm

    outDir = "output/h_za/combine/"
    plot_dir = "/home/submit/jaeyserm/public_html/fccee/h_za/combine_smoothing/"

    if ecm == 240:
        sigs_zh = ['wzp6_ee_nunuH_HZa_ecm240', 'wzp6_ee_eeH_HZa_ecm240', 'wzp6_ee_tautauH_HZa_ecm240', 'wzp6_ee_ccH_HZa_ecm240', 'wzp6_ee_bbH_HZa_ecm240', 'wzp6_ee_qqH_HZa_ecm240', 'wzp6_ee_ssH_HZa_ecm240', 'wzp6_ee_mumuH_HZa_ecm240']
        bkgs = ['wz3p6_ee_gammagamma_ecm240', 'wz3p6_ee_uu_ecm240', 'wz3p6_ee_dd_ecm240', 'wz3p6_ee_cc_ecm240', 'wz3p6_ee_ss_ecm240', 'wz3p6_ee_bb_ecm240', 'wz3p6_ee_tautau_ecm240', 'wz3p6_ee_mumu_ecm240', 'wz3p6_ee_nunu_ecm240', 'wz3p6_ee_ee_Mee_30_150_ecm240', 'p8_ee_ZZ_ecm240']

    if ecm == 365:
        sigs_zh = ['wzp6_ee_numunumuH_HZa_ecm365', 'wzp6_ee_eeH_HZa_ecm365', 'wzp6_ee_tautauH_HZa_ecm365', 'wzp6_ee_ccH_HZa_ecm365', 'wzp6_ee_bbH_HZa_ecm365', 'wzp6_ee_qqH_HZa_ecm365', 'wzp6_ee_ssH_HZa_ecm365', 'wzp6_ee_mumuH_HZa_ecm365']
        sigs_vbf = ['wzp6_ee_nuenueVBFH_HZa_ecm365']
        bkgs = ['wz3p6_ee_gammagamma_ecm365', 'wz3p6_ee_uu_ecm365', 'wz3p6_ee_dd_ecm365', 'wz3p6_ee_cc_ecm365', 'wz3p6_ee_ss_ecm365', 'wz3p6_ee_bb_ecm365', 'wz3p6_ee_tautau_ecm365', 'wz3p6_ee_mumu_ecm365', 'wz3p6_ee_nunu_ecm365', 'wz3p6_ee_ee_Mee_30_150_ecm365', 'p8_ee_ZZ_ecm365']



    cats = ['qqvv']
    if ecm == 365:
        cats.append('vbf')
    for cat in cats:
        print("**************************", cat)
        if ecm == 240:
            sigs = sigs_zh
            inputDir = f"output/h_za/histmaker/ecm{ecm}_qqvv/"
            hName, rebin = 'mva_score', 10

            h_sig = getHist(hName, sigs, rebin=rebin)
            h_bkg = getHist(hName, bkgs, rebin=rebin)

            #h_sig = smooth_histogram_gaussian_kernel(h_sig, sigma, f"{cat}_sig")
            #h_bkg = smooth_histogram_gaussian_kernel(h_bkg, sigma, f"{cat}_bkg")

            h_sig.SetName(f"{cat}_sig_zh")
            h_bkg.SetName(f"{cat}_bkg")

            fOut = ROOT.TFile(f"{outDir}/datacard_{cat}.root", "RECREATE")
            h_sig.Write()
            h_bkg.Write()
            h_data = h_sig.Clone(f"{cat}_data")
            h_data.Add(h_bkg)

            sign = h_sig.Integral() / (h_data.Integral())**0.5
            print(f"SIG: {h_sig.Integral()}")
            print(f"DATA: {h_data.Integral()}")
            print(f"SIGNIFICANCE: {sign}")

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
            dc += "process      sig_zh      bkg\n"
            dc += "process      0           1\n"
            dc += "rate         -1          -1\n"
            dc += "####################\n"
            #dc += "dummy lnN    -     1.00001\n"
            dc += "bkg lnN      -           1.01\n"
            #dc += "* autoMCStats 0 0 1\n"


        if ecm == 365:
            inputDir = f"output/h_za/histmaker/ecm{ecm}_{cat}/"
            hName, rebin = 'mva_score', 10
            h_sig_zh = getHist(hName, sigs_zh, rebin=rebin)
            h_sig_vbf = getHist(hName, sigs_vbf, rebin=rebin)
            h_bkg = getHist(hName, bkgs, rebin=rebin)

            h_sig_vbf = removeNegativeBins(h_sig_vbf)


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
            cmd = f"singularity exec --bind /work:/work /work/submit/jaeyserm/software/docker/cmssw_cc7.sif bash -c 'PYTHONPATH=''; dir=$(pwd); cd /work/submit/jaeyserm/wmass/CMSSW_10_6_19_patch2/src; source /cvmfs/cms.cern.ch/cmsset_default.sh; cmsenv; cd $dir/{outDir}; text2hdf5.py --X-allow-no-background datacard_{cat}.txt -o ws_{cat}.hdf5; combinetf.py ws_{cat}.hdf5 -o fit_output_{cat}.root -t -1  --expectSignal=1  '"
            os.system(cmd)


    if args.run:
        cards = ' '.join([f'datacard_{x}.txt' for x in cats])
        cmd = f"singularity exec --bind /work:/work /work/submit/jaeyserm/software/docker/cmssw_cc7.sif bash -c 'PYTHONPATH=''; dir=$(pwd); cd /work/submit/jaeyserm/wmass/CMSSW_10_6_19_patch2/src; source /cvmfs/cms.cern.ch/cmsset_default.sh; cmsenv; cd $dir/{outDir}; combineCards.py {cards} > datacard_combined.txt'"
        os.system(cmd)

        cmd = f"singularity exec --bind /work:/work /work/submit/jaeyserm/software/docker/cmssw_cc7.sif bash -c 'PYTHONPATH=''; dir=$(pwd); cd /work/submit/jaeyserm/wmass/CMSSW_10_6_19_patch2/src; source /cvmfs/cms.cern.ch/cmsset_default.sh; cmsenv; cd $dir/{outDir}; text2hdf5.py --X-allow-no-background datacard_combined.txt -o ws_combined.hdf5; combinetf.py ws_combined.hdf5 -o fit_output_combined.root -t -1  --expectSignal=1 '" # 
        os.system(cmd)

        for cat in cats+['combined']:
            fIn = ROOT.TFile(f"{outDir}/fit_output_{cat}.root")
            tree = fIn.Get("fitresults")
            tree.GetEntry(0)
            status = tree.status
            errstatus = tree.errstatus

            if ecm == 365:
                mu_zh = tree.sig_zh_mu
                mu_zh_err = tree.sig_zh_mu_err*100.
                mu_vbf = tree.sig_vbf_mu
                mu_vbf_err = tree.sig_vbf_mu_err*100.
                print(f"{cat}\tzh={mu_zh_err:.3f}\tvbf={mu_vbf_err:.3f} {status}")
            else:
                mu_zh = tree.sig_zh_mu
                mu_zh_err = tree.sig_zh_mu_err*100.
                print(f"{cat}\tzh={mu_zh_err:.3f} {status}")



