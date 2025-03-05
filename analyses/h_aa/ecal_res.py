import sys
import os
import math
import copy
import array

import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

sys.path.insert(0, f'{os.path.dirname(os.path.realpath(__file__))}/../../python')
import plotter


if __name__ == "__main__":

    ecm, lumi = 240, "10.8"
    ecm, lumi = 365, "3"

    leg = ROOT.TLegend(0.2, 0.7, 0.9, 0.9)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.035)
    leg.SetMargin(0.2)

    '''
    1   3.424   12.563  14.373
    2   3.507   12.694  14.514
    5   3.967   13.553  15.272
    10  5.116   15.942  16.950
    25  9.868   25.536  22.051
    50  15.750  33.689  31.972
    '''

    x_values = array.array('d', [1, 2, 5, 10, 25])
    y_values_zh_240 = array.array('d', [3.424, 3.5074, 3.967, 5.116, 9.868])
    y_values_zh_365 = array.array('d', [12.563, 12.694, 13.553, 15.942, 25.536])
    y_values_vbf_365 = array.array('d', [14.373, 14.514, 15.272, 18.0, 31.972]) # change 16.950 to 18, as 16.95 interferes a bit with the backgrounds



    graph_zh_240 = ROOT.TGraph(len(x_values), x_values, y_values_zh_240)
    graph_zh_240.SetMarkerStyle(20)
    graph_zh_240.SetMarkerSize(1)
    graph_zh_240.SetLineColor(ROOT.kBlack)
    graph_zh_240.SetMarkerColor(ROOT.kBlack)
    graph_zh_240.SetLineWidth(2)
    leg.AddEntry(graph_zh_240, "ZH #sqrt{s}=240 GeV, 10.8 ab^{#minus1}", "LP")

    graph_zh_365 = ROOT.TGraph(len(x_values), x_values, y_values_zh_365)
    graph_zh_365.SetMarkerStyle(20)
    graph_zh_365.SetMarkerSize(1)
    graph_zh_365.SetLineColor(ROOT.kRed)
    graph_zh_365.SetMarkerColor(ROOT.kRed)
    graph_zh_365.SetLineWidth(2)
    leg.AddEntry(graph_zh_365, "ZH #sqrt{s}=365 GeV, 3 ab^{#minus1}", "LP")

    graph_vbf_365 = ROOT.TGraph(len(x_values), x_values, y_values_vbf_365)
    graph_vbf_365.SetMarkerStyle(20)
    graph_vbf_365.SetMarkerSize(1)
    graph_vbf_365.SetLineColor(ROOT.kBlue)
    graph_vbf_365.SetMarkerColor(ROOT.kBlue)
    graph_vbf_365.SetLineWidth(2)
    leg.AddEntry(graph_vbf_365, "VBF #sqrt{s}=365 GeV, 3 ab^{#minus1}", "LP")



    cfg = {
        'logy'              : False,
        'logx'              : False,

        'xmin'              : 0,
        'xmax'              : 30,
        'ymin'              : 0,
        'ymax'              : 40,
            
        'xtitle'            : "IDEA ECAL stochastic term (%)",
        'ytitle'            : "Uncertainty (%)",
            
        'topRight'          : f"",
        'topLeft'           : "#bf{FCC-ee} #scale[0.7]{#it{Simulation}}",

    }

    plotter.cfg = cfg
    canvas = plotter.canvas()
    dummy = plotter.dummy()
    dummy.Draw("HIST") 
    graph_zh_240.Draw("SAME LP")
    graph_zh_365.Draw("SAME LP")
    graph_vbf_365.Draw("SAME LP")

    leg.Draw("SAME")
    canvas.SetGrid()
    canvas.Modify()
    canvas.Update()

    plotter.aux()
    ROOT.gPad.SetTicks()
    ROOT.gPad.RedrawAxis()

    canvas.SaveAs("/home/submit/jaeyserm/public_html/fccee/h_aa/ecal_variations.png")
    canvas.SaveAs("/home/submit/jaeyserm/public_html/fccee/h_aa/ecal_variations.pdf")
    canvas.Close()