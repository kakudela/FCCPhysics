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



def comparison(hName, outName, procs, xMin=0, xMax=100, yMin=1, yMax=1e5, xLabel="xlabel", yLabel="Events", logX=False, logY=True, rebin=1, legPos=[0.3, 0.75, 0.9, 0.9]):
    st = ROOT.THStack()
    st.SetName("stack")
    
    colors = [ROOT.kBlue, ROOT.kRed]

    leg = ROOT.TLegend(legPos[0], legPos[1], legPos[2], legPos[3])
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.03)
    leg.SetMargin(0.2)

    hists = []
    for i,proc in enumerate(procs):
        hist = getHist(hName, [proc], rebin)
        hist.SetName(proc)
        hist.SetLineColor(colors[i])
        hist.SetLineWidth(2)
        hist.SetLineStyle(1)
        leg.AddEntry(hist, proc, "L")
        
        hist.Scale(1./hist.Integral())
        hists.append(hist)

    if yMax < 0:
        if logY:
            yMax = math.ceil(max([h.GetMaximum() for h in hists])*10000)/10.
        else:
            yMax = 1.2*max([h.GetMaximum() for h in hists])

    cfg = {

        'logy'              : logY,
        'logx'              : logX,
        
        'xmin'              : xMin,
        'xmax'              : xMax,
        'ymin'              : yMin,
        'ymax'              : yMax,
            
        'xtitle'            : xLabel,
        'ytitle'            : yLabel,
            
        'topRight'          : f"#sqrt{{s}} = {ecm} GeV, {lumi} ab^{{#minus1}}",
        'topLeft'           : "#bf{FCC-ee} #scale[0.7]{#it{Simulation}}",

    }

    plotter.cfg = cfg
    canvas = plotter.canvas()
    dummy = plotter.dummy()
    dummy.Draw("HIST") 

    for hist in hists:
        hist.Draw("SAME HIST")

    leg.Draw("SAME")
    
    canvas.SetGrid()
    canvas.Modify()
    canvas.Update()

    plotter.aux()
    ROOT.gPad.SetTicks()
    ROOT.gPad.RedrawAxis()

    canvas.SaveAs(f"{outDir}/{outName}.png")
    canvas.SaveAs(f"{outDir}/{outName}.pdf")
    canvas.Close()




if __name__ == "__main__":

    inputDir = "output/kkmcee_whizard/histmaker/"
    outDir = "/home/submit/jaeyserm/public_html/fccee/kkmcee_whizard/plots/"
    
    ecm, lumi = 240, "10.8"


    procs = ['kkmcee_ee_mumu_ecm240', 'wzp6_ee_mumu_ecm240']


    comparison("gen_photons_no", "gen_photons_no", procs, xMin=0, xMax=10, yMax=-1, xLabel="Photon multiplicity", yLabel="Events", logY=False)
    comparison("gen_photons_p", "gen_photons_p", procs, xMin=0, xMax=0.1, yMax=-1, xLabel="Photon momentum (GeV)", yLabel="Events", logY=False)
    comparison("gen_photons_theta", "gen_photons_theta", procs, xMin=0, xMax=3.15, yMax=-1, xLabel="Photon theta (rad)", yLabel="Events", logY=False)

    comparison("m_ll", "m_ll", procs, xMin=0, xMax=250, yMax=-1, xLabel="m_ll", yLabel="Events", logY=False)
    comparison("p_ll", "p_ll", procs, xMin=0, xMax=250, yMax=-1, xLabel="p_ll", yLabel="Events", logY=False)

    comparison("leptons_p", "leptons_p", procs, xMin=0, xMax=250, yMax=-1, xLabel="Leptons momentum (GeV)", yLabel="Events", logY=False)
    comparison("visibleEnergy", "visibleEnergy", procs, xMin=0, xMax=250, yMax=-1, xLabel="visibleEnergy", yLabel="Events", logY=False)


    comparison("photon1_p", "photon1_p", procs, xMin=0, xMax=150, yMin=0.0001, yMax=-1, xLabel="photon1_p", yLabel="Events", logY=True)
    comparison("photon2_p", "photon2_p", procs, xMin=0, xMax=150, yMin=0.0001, yMax=-1, xLabel="photon2_p", yLabel="Events", logY=True)



    comparison("photon1_p_rr", "photon1_p_rr", procs, xMin=0, xMax=150, yMin=0.0001, yMax=-1, xLabel="photon1_p_rr", yLabel="Events", logY=True)
    comparison("photon2_p_rr", "photon2_p_rr", procs, xMin=0, xMax=150, yMin=0.0001, yMax=-1, xLabel="photon2_p_rr", yLabel="Events", logY=True)
    comparison("p_ll_rr", "p_ll_rr", procs, xMin=0, xMax=250, yMax=-1, xLabel="p_ll_rr", yLabel="Events", logY=False)
    comparison("visibleEnergy_rr", "visibleEnergy_rr", procs, xMin=0, xMax=250, yMax=-1, xLabel="visibleEnergy_rr", yLabel="Events", logY=False)
