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
    #hist.Scale(3./10.8 * 1.05/1.39)
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

    ecm = 365
    cat = "vbf"

    inputDir = f"output/h_za/histmaker/ecm{ecm}_{cat}//"
    outDir = "output/h_za/combine/"
    plot_dir = "/home/submit/jaeyserm/public_html/fccee/h_za/combine_smoothing/"
    rebin = 1 # 100 MVA bins
    sigma = 1

    if ecm == 240:
        sigs = ['wzp6_ee_nunuH_HZa_ecm240', 'wzp6_ee_eeH_HZa_ecm240', 'wzp6_ee_tautauH_HZa_ecm240', 'wzp6_ee_ccH_HZa_ecm240', 'wzp6_ee_bbH_HZa_ecm240', 'wzp6_ee_qqH_HZa_ecm240', 'wzp6_ee_ssH_HZa_ecm240', 'wzp6_ee_mumuH_HZa_ecm240']
        bkgs = ['wz3p6_ee_gammagamma_ecm240', 'wz3p6_ee_uu_ecm240', 'wz3p6_ee_dd_ecm240', 'wz3p6_ee_cc_ecm240', 'wz3p6_ee_ss_ecm240', 'wz3p6_ee_bb_ecm240', 'wz3p6_ee_tautau_ecm240', 'wz3p6_ee_mumu_ecm240', 'wz3p6_ee_nunu_ecm240', 'wz3p6_ee_ee_Mee_30_150_ecm240', 'p8_ee_ZZ_ecm240']

    if ecm == 365:
        sigs = ['wzp6_ee_nunuH_HZa_ecm240', 'wzp6_ee_eeH_HZa_ecm240', 'wzp6_ee_tautauH_HZa_ecm240', 'wzp6_ee_ccH_HZa_ecm240', 'wzp6_ee_bbH_HZa_ecm240', 'wzp6_ee_qqH_HZa_ecm240', 'wzp6_ee_ssH_HZa_ecm240', 'wzp6_ee_mumuH_HZa_ecm240']
        sigs = ['wzp6_ee_nuenueVBFH_HZa_ecm365']
        bkgs = ['wz3p6_ee_gammagamma_ecm365', 'wz3p6_ee_uu_ecm365', 'wz3p6_ee_dd_ecm365', 'wz3p6_ee_cc_ecm365', 'wz3p6_ee_ss_ecm365', 'wz3p6_ee_bb_ecm365', 'wz3p6_ee_tautau_ecm365', 'wz3p6_ee_mumu_ecm365', 'wz3p6_ee_nunu_ecm365', 'wz3p6_ee_ee_Mee_30_150_ecm365', 'p8_ee_ZZ_ecm365']

    #bkgs = ['p8_ee_ZZ_ecm240']

    hists = ["hqqa_final_1D_mva", "hvva_final_1D_mva"] # 28.976/34.240/23.235
    hists = ["hqqa_final_chi2", "hvva_final_chi2"] # 28.733/34.454/23.147
    hists = ['qqvv_H_m'] # chi2 pairing
    ## seems no difference between chi2 and MVA?
    hists = ["hvva_mva_score", "hqqa_mva_score"] 
    #hists = ["hqqa_final_2D_mva", "hvva_final_2D_mva"]
    
    #hists = ["hqqa_final_chi2", "hvva_final_chi2"] # 23%
    #hists = ["hqqa_final_1D_mva", "hvva_final_1D_mva"] # 23.2 %
    #hists = ["qqvv_mva_score"]

    # qqqq tests
    #hists = ["qqqq_H_m_final"] # 55.331
    #hists = ["qqqq_H_m_hz1", "qqqq_H_m_hz2"] # 69.121/73.117/55.855
    
    
    #hists = ["hqqa_mass_difference_qqa_vv", "hvva_mass_difference_vva_qq"] # 38.942/43.131/32.034
    #hists = ["qqqq_mass_difference_Z1a_Z2", "qqqq_mass_difference_Z2a_Z1"] # 67.684/69.297/56.533
    
    # pairing_chi2
    hists = ["qqvv_H_m"] # 23.242
    hists = ["qqvv_mass_difference_H_Z"] # 34.509
    hists = ["qqvv_mva_score_split_chi2"] # 12.134 /// 11.481 (full ZZ stats) --> factor of 2 with MVA
    
    #hists = ["qqvv_hqqa_mva_score", "qqvv_hvva_mva_score"]

    # splitting_chi2
    #hists = ["qqvv_hqqa_final_chi2", "qqvv_hvva_final_chi2"] # 27.225/34.781/22.598
    hists = ["qqvv_qqa_m_vva_m"]

    # splitting_mva
    #hists = ["qqvv_hqqa_final_1D_mva", "qqvv_hvva_final_1D_mva"] # 27.254/34.996/22.720
    
    if ecm == 365 and cat == "vbf":
        hists = ['qqa_m']

    for hName in hists:
        h_sig = getHist(hName, sigs, rebin=rebin)
        h_bkg = getHist(hName, bkgs, rebin=rebin)

        #h_sig = unroll(h_sig)
        #h_bkg = unroll(h_bkg)
        #h_sig = smooth_histogram_gaussian_kernel(h_sig, sigma, f"{hName}_sig")
        #h_bkg = smooth_histogram_gaussian_kernel(h_bkg, sigma, f"{hName}_bkg")

        h_sig.SetName(f"{hName}_sig")
        h_bkg.SetName(f"{hName}_bkg")

        fOut = ROOT.TFile(f"{outDir}/datacard_{hName}.root", "RECREATE")
        h_sig.Write()
        h_bkg.Write()
        h_data = h_sig.Clone(f"{hName}_data")
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
        #dc += "dummy lnN    -     1.00001\n"
        dc += "bkg lnN      -           1.01\n"
        #dc += "* autoMCStats 0 0 1\n"


        f = open(f"{outDir}/datacard_{hName}.txt", 'w')
        f.write(dc)
        f.close()

        print(dc)

        if args.run:
            cmd = f"singularity exec --bind /work:/work /work/submit/jaeyserm/software/docker/cmssw_cc7.sif bash -c 'PYTHONPATH=''; dir=$(pwd); cd /work/submit/jaeyserm/wmass/CMSSW_10_6_19_patch2/src; source /cvmfs/cms.cern.ch/cmsset_default.sh; cmsenv; cd $dir/{outDir}; text2hdf5.py --X-allow-no-background datacard_{hName}.txt -o ws_{hName}.hdf5; combinetf.py ws_{hName}.hdf5 -o fit_output_{hName}.root -t 0  --expectSignal=1'"
            os.system(cmd)
    #quit()

    if args.run:
        cards = ' '.join([f'datacard_{x}.txt' for x in hists])
        cmd = f"singularity exec --bind /work:/work /work/submit/jaeyserm/software/docker/cmssw_cc7.sif bash -c 'PYTHONPATH=''; dir=$(pwd); cd /work/submit/jaeyserm/wmass/CMSSW_10_6_19_patch2/src; source /cvmfs/cms.cern.ch/cmsset_default.sh; cmsenv; cd $dir/{outDir}; combineCards.py {cards} > datacard_combined.txt'"
        os.system(cmd)

        cmd = f"singularity exec --bind /work:/work /work/submit/jaeyserm/software/docker/cmssw_cc7.sif bash -c 'PYTHONPATH=''; dir=$(pwd); cd /work/submit/jaeyserm/wmass/CMSSW_10_6_19_patch2/src; source /cvmfs/cms.cern.ch/cmsset_default.sh; cmsenv; cd $dir/{outDir}; text2hdf5.py --X-allow-no-background datacard_combined.txt -o ws_combined.hdf5; combinetf.py ws_combined.hdf5 -o fit_output_combined.root -t 0  --expectSignal=1 '" # --binByBinStat
        os.system(cmd)

        for hist in hists+['combined']:
            fIn = ROOT.TFile(f"{outDir}/fit_output_{hist}.root")
            tree = fIn.Get("fitresults")
            tree.GetEntry(0)
            status = tree.status
            errstatus = tree.errstatus

            mu_zh = tree.sig_mu
            mu_zh_err = tree.sig_mu_err*100.
            print(f"{hist}\tzh={mu_zh_err:.3f} {status}")


'''
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

'''