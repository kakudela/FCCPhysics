
import sys,array,ROOT,math,os,copy
import numpy as np

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)


def compute_res(input_file, hist_name, output_name, plotGauss=True):


    fIn = ROOT.TFile(input_file)
    hist = fIn.Get(hist_name)

    rebin = 1
    hist = hist.Rebin(rebin)

    probabilities = np.array([0.001, 0.999, 0.84, 0.16], dtype='d')


    # compute quantiles
    quantiles = np.array([0.0, 0.0, 0.0, 0.0], dtype='d')
    hist.GetQuantiles(4, quantiles, probabilities)
    xMin, xMax = min([quantiles[0], -quantiles[1]]), max([-quantiles[0], quantiles[1]])
    res = 100.*0.5*(quantiles[2] - quantiles[3])

    # compute RMS
    rms, rms_err = hist.GetRMS()*100., hist.GetRMSError()*100.

    # fit with Gauss
    gauss = ROOT.TF1("gauss2", "gaus", xMin, xMax)
    gauss.SetParameter(0, hist.Integral())
    gauss.SetParameter(1, hist.GetMean())
    gauss.SetParameter(2, hist.GetRMS())
    hist.Fit("gauss2", "R")

    sigma, sigma_err = gauss.GetParameter(2)*100., gauss.GetParError(2)*100.

    gauss.SetLineColor(ROOT.kRed)
    gauss.SetLineWidth(3)


    yMin, yMax = 0, 1.3*hist.GetMaximum()
    canvas = ROOT.TCanvas("canvas", "", 1000, 1000)
    canvas.SetTopMargin(0.055)
    canvas.SetRightMargin(0.05)
    canvas.SetLeftMargin(0.15)
    canvas.SetBottomMargin(0.11)

    dummy = ROOT.TH1D("h", "h", 1, xMin, xMax)
    dummy.GetXaxis().SetTitle("(p_{reco} #minus p_{gen})/p_{gen}")
    dummy.GetXaxis().SetRangeUser(xMin, xMax)

    dummy.GetXaxis().SetTitleFont(43)
    dummy.GetXaxis().SetTitleSize(40)
    dummy.GetXaxis().SetLabelFont(43)
    dummy.GetXaxis().SetLabelSize(35)

    dummy.GetXaxis().SetTitleOffset(1.2*dummy.GetXaxis().GetTitleOffset())
    dummy.GetXaxis().SetLabelOffset(1.2*dummy.GetXaxis().GetLabelOffset())

    dummy.GetYaxis().SetTitle("Events / bin")
    dummy.GetYaxis().SetRangeUser(yMin, yMax)
    dummy.SetMaximum(yMax)
    dummy.SetMinimum(yMin)

    dummy.GetYaxis().SetTitleFont(43)
    dummy.GetYaxis().SetTitleSize(40)
    dummy.GetYaxis().SetLabelFont(43)
    dummy.GetYaxis().SetLabelSize(35)

    dummy.GetYaxis().SetTitleOffset(1.7*dummy.GetYaxis().GetTitleOffset())
    dummy.GetYaxis().SetLabelOffset(1.4*dummy.GetYaxis().GetLabelOffset())

    dummy.Draw("HIST")
    hist.Draw("SAME HIST")
    if plotGauss:
        gauss.Draw("SAME")

    canvas.SetGrid()
    ROOT.gPad.SetTickx()
    ROOT.gPad.SetTicky()
    ROOT.gPad.RedrawAxis()

    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.035)
    latex.SetTextColor(1)
    latex.SetTextFont(42)
    latex.DrawLatex(0.2, 0.9, f"Mean/RMS(#times 100) = {hist.GetMean():.4f}/{rms:.4f}")
    latex.DrawLatex(0.2, 0.85, f"Resolution = {res:.4f} %")
    if plotGauss:
        latex.DrawLatex(0.2, 0.80, f"Gauss #mu/#sigma(#times 100) = {gauss.GetParameter(1):.4f}/{sigma:.4f}")

    canvas.SaveAs(f"{output_name}.png")
    canvas.SaveAs(f"{output_name}.pdf")
    canvas.Close()

    del gauss
    return rms, rms_err, sigma, sigma_err, res



if __name__ == "__main__":

    input_file, output_name = "output/IDEA_2T_Zmumu_ecm240.root", "resolution_IDEA_2T_Zmumu_ecm240"
    input_file, output_name = "output/IDEA_3T_Zmumu_ecm240.root", "resolution_IDEA_3T_Zmumu_ecm240"
    #input_file, output_name = "output/CLD_2T_Zmumu_ecm240.root", "resolution_CLD_2T_Zmumu_ecm240"
    hist_name = "muons_res_p"


    compute_res(input_file, hist_name, output_name)
