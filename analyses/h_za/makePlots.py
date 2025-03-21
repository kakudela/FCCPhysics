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


import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--cat", type=str, help="Category (qqvv, qqqq)", default="qqvv")
parser.add_argument("--ecm", type=int, help="Center-of-mass energy", default=240)
args = parser.parse_args()






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

def makeCutFlow(hName="cutFlow", cuts=[], labels=[], sig_scale=1.0, yMin=1e6, yMax=1e10):

    leg = ROOT.TLegend(.55, 0.99-(len(procs))*0.06, .99, .90)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.03)
    leg.SetMargin(0.2)

    hists_yields = []
    significances = []
    h_sig = getHist(hName, procs_cfg[procs[0]])

    hists_yields.append(copy.deepcopy(h_sig))
    h_sig.Scale(sig_scale)
    h_sig.SetLineColor(procs_colors[procs[0]])
    h_sig.SetLineWidth(4)
    h_sig.SetLineStyle(1)
    if sig_scale != 1:
        leg.AddEntry(h_sig, f"{procs_labels[procs[0]]} (#times {int(sig_scale)})", "L")
    else:
        leg.AddEntry(h_sig, procs_labels[procs[0]], "L")

    st = ROOT.THStack()
    st.SetName("stack")
    h_bkg_tot = None
    for i,bkg in enumerate(procs[1:]):
        h_bkg = getHist(hName, procs_cfg[bkg])

        if h_bkg_tot == None: h_bkg_tot = h_bkg.Clone("h_bkg_tot")
        else: h_bkg_tot.Add(h_bkg)
        
        h_bkg.SetFillColor(procs_colors[bkg])
        h_bkg.SetLineColor(ROOT.kBlack)
        h_bkg.SetLineWidth(1)
        h_bkg.SetLineStyle(1)

        leg.AddEntry(h_bkg, procs_labels[bkg], "F")
        st.Add(h_bkg)
        hists_yields.append(h_bkg)

    h_bkg_tot.SetLineColor(ROOT.kBlack)
    h_bkg_tot.SetLineWidth(2)

    for i,cut in enumerate(cuts):
        nsig = h_sig.GetBinContent(i+1) / sig_scale ## undo scaling
        nbkg = 0
        for j,histProc in enumerate(hists_yields):
            nbkg = nbkg + histProc.GetBinContent(i+1)
            print(histProc.GetBinContent(i+1))
        if (nsig+nbkg) == 0:
            print(f"Cut {cut} zero yield sig+bkg")
            s = -1
        else:
            s = nsig / (nsig + nbkg)**0.5
        print(i, cut, s)
        significances.append(s)

    ########### PLOTTING ###########
    cfg = {
        'logy'              : True,
        'logx'              : False,

        'xmin'              : 0,
        'xmax'              : len(cuts),
        'ymin'              : yMin,
        'ymax'              : yMax ,

        'xtitle'            : "",
        'ytitle'            : "Events",

        'topRight'          : f"#sqrt{{s}} = {ecm} GeV, {lumi} ab^{{#minus1}}",
        'topLeft'           : "#bf{FCC-ee} #scale[0.7]{#it{Simulation}}",
        }

    plotter.cfg = cfg

    canvas = plotter.canvas()
    canvas.SetGrid()
    canvas.SetTicks()
    dummy = plotter.dummy(len(cuts))
    dummy.GetXaxis().SetLabelSize(0.75*dummy.GetXaxis().GetLabelSize())
    dummy.GetXaxis().SetLabelOffset(1.3*dummy.GetXaxis().GetLabelOffset())
    for i,label in enumerate(labels): dummy.GetXaxis().SetBinLabel(i+1, label)
    dummy.GetXaxis().LabelsOption("u")
    dummy.Draw("HIST")

    st.Draw("SAME HIST")
    h_bkg_tot.Draw("SAME HIST")
    h_sig.Draw("SAME HIST")
    leg.Draw("SAME")

    plotter.aux()
    canvas.RedrawAxis()
    canvas.Modify()
    canvas.Update()
    canvas.Draw()
    canvas.SaveAs(f"{outDir}/{hName}.png")
    canvas.SaveAs(f"{outDir}/{hName}.pdf")

    out_orig = sys.stdout
    with open(f"{outDir}{hName}.txt", 'w') as f:
        sys.stdout = f

        formatted_row = '{:<10} {:<25} ' + ' '.join(['{:<25}']*len(procs))
        print(formatted_row.format(*(["Cut", "Significance"]+procs)))
        print(formatted_row.format(*(["----------"]+["-----------------------"]*(len(procs)+1))))
        for i,cut in enumerate(cuts):
            row = ["Cut %d"%i, "%.3f"%significances[i]]
            for j,histProc in enumerate(hists_yields):
                yield_, err = histProc.GetBinContent(i+1), histProc.GetBinError(i+1)
                row.append("%.4e +/- %.2e" % (yield_, err))

            print(formatted_row.format(*row))
    sys.stdout = out_orig

def makeCutFlowOld():

    totEntries = 1 + len(bkgs)
    #leg = ROOT.TLegend(.5, 1.0-totEntries*0.06, .92, .90)
    leg = ROOT.TLegend(.45, 0.99-(len(bkgs)+2)*0.055, .95, .90)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.03)
    leg.SetMargin(0.2)

    hists_yields = []
    h_sig = getHist("cutFlow", sigs)

    hists_yields.append(copy.deepcopy(h_sig))
    h_sig.Scale(sig_scale)
    h_sig.SetLineColor(ROOT.TColor.GetColor("#BF2229"))
    h_sig.SetLineWidth(4)
    h_sig.SetLineStyle(1)
    leg.AddEntry(h_sig, sig_legend, "L")

    # Get all bkg histograms
    st = ROOT.THStack()
    st.SetName("stack")
    h_bkg_tot = None
    for i,bkg in enumerate(bkgs):
        h_bkg = getHist("cutFlow", bgks_cfg[bkg])

        if h_bkg_tot == None: h_bkg_tot = h_bkg.Clone("h_bkg_tot")
        else: h_bkg_tot.Add(h_bkg)
        
        h_bkg.SetFillColor(bkgs_colors[i])
        h_bkg.SetLineColor(ROOT.kBlack)
        h_bkg.SetLineWidth(1)
        h_bkg.SetLineStyle(1)

        leg.AddEntry(h_bkg, bkgs_legends[i], "F")
        st.Add(h_bkg)
        hists_yields.append(h_bkg)

    h_bkg_tot.SetLineColor(ROOT.kBlack)
    h_bkg_tot.SetLineWidth(2)


    ########### PLOTTING ###########
    cfg = {

        'logy'              : True,
        'logx'              : False,
        
        'xmin'              : 0,
        'xmax'              : len(cuts),
        'ymin'              : 1e2,
        'ymax'              : 1e10 ,
            
        'xtitle'            : "",
        'ytitle'            : "Events",
            
        'topRight'          : "#sqrt{s} = 240 GeV, 7.2 ab^{#minus1}", 
        'topLeft'           : "#bf{FCC-ee} #scale[0.7]{#it{Simulation}}",
        }

    plotter.cfg = cfg

    canvas = plotter.canvas()
    canvas.SetGrid()
    canvas.SetTicks()
    dummy = plotter.dummy(len(cuts))
    dummy.GetXaxis().SetLabelSize(0.8*dummy.GetXaxis().GetLabelSize())
    dummy.GetXaxis().SetLabelOffset(1.3*dummy.GetXaxis().GetLabelOffset())
    for i,label in enumerate(labels): dummy.GetXaxis().SetBinLabel(i+1, label)
    dummy.GetXaxis().LabelsOption("u")
    dummy.Draw("HIST")

    st.Draw("SAME HIST")
    h_bkg_tot.Draw("SAME HIST")
    h_sig.Draw("SAME HIST")
    leg.Draw("SAME")

    plotter.aux()
    canvas.RedrawAxis()
    canvas.Modify()
    canvas.Update()
    canvas.Draw()
    canvas.SaveAs(f"{outDir}/cutFlow.png")
    canvas.SaveAs(f"{outDir}/cutFlow.pdf")
    
    out_orig = sys.stdout
    with open(f"{outDir}cutFlow.txt", 'w') as f:
        sys.stdout = f

        formatted_row = '{:<10} {:<25} {:<25} {:<25}' # adapt to #bkgs
        print(formatted_row.format(*(["Cut", "Signal"]+bkgs)))
        print(formatted_row.format(*(["----------"]+["-----------------------"]*5)))
        for i,cut in enumerate(cuts):
            row = ["Cut %d"%i]
            for j,histProc in enumerate(hists_yields):
                yield_, err = histProc.GetBinContent(i+1), histProc.GetBinError(i+1)
                row.append("%.2e +/- %.2e" % (yield_, err))

            print(formatted_row.format(*row))
    sys.stdout = out_orig

def comparison(hName, outName, xMin=0, xMax=100, yMin=1, yMax=1e5, xLabel="xlabel", yLabel="Events", logX=False, logY=True, rebin=1, legPos=[0.3, 0.75, 0.9, 0.9]):
    st = ROOT.THStack()
    st.SetName("stack")

    leg = ROOT.TLegend(legPos[0], legPos[1], legPos[2], legPos[3])
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.03)
    leg.SetMargin(0.2)

    p8 = getHist(hName, ['p8_ee_Zqq_ecm240'], rebin)
    p8.SetName("p8")
    p8.SetLineColor(ROOT.kRed)
    p8.SetLineWidth(2)
    p8.SetLineStyle(1)
    
    leg.AddEntry(p8, "Pythia8", "L")

    kkmc = getHist(hName, ['kkmcee_ee_uu_ecm240', 'kkmcee_ee_dd_ecm240', 'kkmcee_ee_cc_ecm240', 'kkmcee_ee_ss_ecm240', 'kkmcee_ee_bb_ecm240'], rebin)
    kkmc.SetName("kkmc")
    kkmc.SetLineColor(ROOT.kBlue)
    kkmc.SetLineWidth(2)
    kkmc.SetLineStyle(1)
    
    leg.AddEntry(kkmc, "KKMC", "L")
    
    if False:
        kkmc.Scale(1./kkmc.Integral())
        p8.Scale(1./p8.Integral())

    if yMax < 0:
        if logY:
            yMax = math.ceil(max([p8.GetMaximum(), kkmc.GetMaximum()])*10000)/10.
        else:
            yMax = 1.2*max([p8.GetMaximum(), kkmc.GetMaximum()])

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

    kkmc.Draw("SAME HIST")
    p8.Draw("SAME HIST")

    leg.Draw("SAME")
    
    canvas.SetGrid()
    canvas.Modify()
    canvas.Update()

    plotter.aux()
    ROOT.gPad.SetTicks()
    ROOT.gPad.RedrawAxis()

    canvas.SaveAs(f"{outDir_comparison}/{outName}.png")
    canvas.SaveAs(f"{outDir_comparison}/{outName}.pdf")
    canvas.Close()



def makePlotOld(hName, outName, xMin=0, xMax=100, yMin=1, yMax=1e5, xLabel="xlabel", yLabel="Events", logX=False, logY=True, rebin=1, legPos=[0.3, 0.75, 0.9, 0.9]):


    st = ROOT.THStack()
    st.SetName("stack")

    leg = ROOT.TLegend(legPos[0], legPos[1], legPos[2], legPos[3])
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.03)
    leg.SetMargin(0.2)

    h_sig = getHist(hName, sigs, rebin)
    if "TH2" in h_sig.ClassName(): h_sig = h_sig.ProjectionX("h_sig")
    h_sig.SetLineColor(ROOT.TColor.GetColor("#BF2229"))
    h_sig.SetLineWidth(3)
    h_sig.SetLineStyle(1)
    h_sig.Scale(sig_scale)
    leg.AddEntry(h_sig, sig_legend, "L")

    st = ROOT.THStack()
    st.SetName("stack")
    h_bkg_tot = None
    for i,bkg in enumerate(bkgs):

        hist = getHist(hName, bgks_cfg[bkg], rebin)
        if "TH2" in hist.ClassName(): hist = hist.ProjectionX()
        hist.SetName(bkg)
        hist.SetFillColor(bkgs_colors[i])
        hist.SetLineColor(ROOT.kBlack)
        hist.SetLineWidth(1)
        hist.SetLineStyle(1)
        #hist.Scale(0.88)
        leg.AddEntry(hist, bkgs_legends[i], "F")
        st.Add(hist)
        if h_bkg_tot == None:
            h_bkg_tot = copy.deepcopy(hist)
            h_bkg_tot.SetName("h_bkg_tot")
        else: h_bkg_tot.Add(hist)

    if yMax < 0:
        if logY:
            yMax = math.ceil(max([h_bkg_tot.GetMaximum(), h_sig.GetMaximum()])*10000)/10.
        else:
            yMax = 1.4*max([h_bkg_tot.GetMaximum(), h_sig.GetMaximum()])

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
    st.Draw("HIST SAME")

    '''
    hTot_err = hTot.Clone("hTot_err")
    hTot_err.SetFillColor(ROOT.kBlack)
    hTot_err.SetMarkerColor(ROOT.kBlack)
    hTot_err.SetFillStyle(3004)
    leg.AddEntry(hTot_err, "Stat. Unc.", "F")
    '''

    h_bkg_tot.SetLineColor(ROOT.kBlack)
    h_bkg_tot.SetFillColor(0)
    h_bkg_tot.SetLineWidth(2)
    #hTot_err.Draw("E2 SAME")
    h_bkg_tot.Draw("HIST SAME")
    h_sig.Draw("HIST SAME")
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



def makePlot(hName, outName="", xMin=0, xMax=100, yMin=1, yMax=1e5, xLabel="xlabel", yLabel="Events", logX=False, logY=True, rebin=1, legPos=[0.3, 0.75, 0.9, 0.9], sig_scale=1, xLabels=[]):

    if outName == "":
        outName = hName

    st = ROOT.THStack()
    st.SetName("stack")

    leg = ROOT.TLegend(.55, 0.99-(len(procs))*0.06, .99, .90)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.03)
    leg.SetMargin(0.2)

    h_sig = getHist(hName, procs_cfg[procs[0]], rebin=rebin)
    h_sig.SetLineColor(procs_colors[procs[0]])
    h_sig.SetLineWidth(3)
    h_sig.SetLineStyle(1)
    h_sig.Scale(sig_scale)
    if sig_scale != 1:
        leg.AddEntry(h_sig, f"{procs_labels[procs[0]]} (#times {int(sig_scale)})", "L")
    else:
        leg.AddEntry(h_sig, procs_labels[procs[0]], "L")

    st = ROOT.THStack()
    st.SetName("stack")
    h_bkg_tot = None
    for i,bkg in enumerate(procs[1:]):
        h_bkg = getHist(hName, procs_cfg[bkg], rebin=rebin)

        if h_bkg_tot == None: h_bkg_tot = h_bkg.Clone("h_bkg_tot")
        else: h_bkg_tot.Add(h_bkg)
        
        h_bkg.SetFillColor(procs_colors[bkg])
        h_bkg.SetLineColor(ROOT.kBlack)
        h_bkg.SetLineWidth(1)
        h_bkg.SetLineStyle(1)

        leg.AddEntry(h_bkg, procs_labels[bkg], "F")
        st.Add(h_bkg)

    if yMax < 0:
        if logY:
            yMax = math.ceil(max([h_bkg_tot.GetMaximum(), h_sig.GetMaximum()])*10000)/10.
        else:
            yMax = 1.4*max([h_bkg_tot.GetMaximum(), h_sig.GetMaximum()])

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
    dummy = plotter.dummy(1 if len(xLabels) == 0 else len(xLabels))
    if len(xLabels) > 0:
        dummy.GetXaxis().SetLabelSize(0.8*dummy.GetXaxis().GetLabelSize())
        dummy.GetXaxis().SetLabelOffset(1.3*dummy.GetXaxis().GetLabelOffset())
        for i,label in enumerate(xLabels): dummy.GetXaxis().SetBinLabel(i+1, label)
        dummy.GetXaxis().LabelsOption("u")
    dummy.Draw("HIST")
    st.Draw("HIST SAME")

    h_bkg_tot.SetLineColor(ROOT.kBlack)
    h_bkg_tot.SetFillColor(0)
    h_bkg_tot.SetLineWidth(2)
    #hTot_err.Draw("E2 SAME")
    h_bkg_tot.Draw("HIST SAME")
    h_sig.Draw("HIST SAME")
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


def significanceold(hName, xMin=-10000, xMax=10000, reverse=False):

    h_sig = getHist(hName, sigs)
    sig_tot = h_sig.Integral()

    bkgs_procs = []
    for bkg in bkgs:
        bkgs_procs.extend(bgks_cfg[bkg])

    h_bkg = getHist(hName, bkgs_procs)
    x, y, l = [], [], []

    for i in range(1, h_sig.GetNbinsX()+1):
        if reverse:
            iStart = 1
            iEnd = i
        else:
            iStart = i
            iEnd = h_sig.GetNbinsX()+1
        center = h_sig.GetBinCenter(i)
        if center > xMax or center < xMin:
            continue
        sig = h_sig.Integral(iStart, iEnd)
        bkg = h_bkg.Integral(iStart, iEnd)
        significance = sig / (sig + bkg)**0.5
        sig_loss = sig / sig_tot
        print(f"{i} {center:.5f} {significance:.5f} {sig_loss:.5f}")
        x.append(center)
        y.append(significance)
        l.append(sig_loss)


    graph_sign = ROOT.TGraph(len(x), array.array('d', x), array.array('d', y))
    graph_l = ROOT.TGraph(len(x), array.array('d', x), array.array('d', l))
    max_y = max(y)
    max_index = y.index(max_y)
    max_x = x[max_index]
    max_l = l[max_index]

    canvas = ROOT.TCanvas("", "", 800, 800)
    graph_sign.SetMarkerStyle(20)
    graph_sign.SetMarkerColor(ROOT.kBlue)
    graph_sign.GetXaxis().SetRangeUser(xMin, xMax)
    graph_sign.Draw("AP")
    canvas.Update()

    # Add a marker for the maximum point
    max_marker = ROOT.TMarker(max_x, max_y, 20)
    max_marker.SetMarkerColor(ROOT.kRed)
    max_marker.SetMarkerSize(1.5)
    max_marker.Draw()



    rightmax = 1.0
    print(ROOT.gPad.GetUymin())
    print(ROOT.gPad.GetUymax())
    scale =  ROOT.gPad.GetUymax()/rightmax
    rightmin = ROOT.gPad.GetUymin()/ROOT.gPad.GetUymax()
    graph_l.Scale(scale)
    graph_l.SetLineColor(ROOT.kRed)
    graph_l.SetLineWidth(2)
    graph_l.Draw("SAME L")


    axis_r = ROOT.TGaxis(ROOT.gPad.GetUxmax(),ROOT.gPad.GetUymin(), ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymax(),rightmin,rightmax,510,"+L")
    axis_r.SetLineColor(ROOT.kRed)
    axis_r.SetLabelColor(ROOT.kRed)
    axis_r.Draw()


    # Add a text box to indicate the maximum value
    text = ROOT.TLatex()
    text.SetTextSize(0.03)
    text.SetTextColor(ROOT.kBlack)
    text.DrawLatexNDC(0.1, 0.92, f"Max: x = {max_x}, y = {max_y:.5f}, signal loss = {max_l:.5f}")


    suffix = "_reverse" if reverse else ""
    canvas.SaveAs(f"{outDir}/significance/{hName}{suffix}.png")
    canvas.SaveAs(f"{outDir}/significance/{hName}{suffix}.pdf")
    canvas.Close()


def significance(hName, xMin=-10000, xMax=10000, reverse=False):

    h_sig = getHist(hName, procs_cfg[procs[0]])
    sig_tot = h_sig.Integral()

    bkgs_procs = []
    for i,bkg in enumerate(procs[1:]):
        bkgs_procs.extend(procs_cfg[bkg])

    h_bkg = getHist(hName, bkgs_procs)
    x, y, l = [], [], []

    for i in range(1, h_sig.GetNbinsX()+1):
        if reverse:
            iStart = 1
            iEnd = i
        else:
            iStart = i
            iEnd = h_sig.GetNbinsX()+1
        center = h_sig.GetBinCenter(i)
        if center > xMax or center < xMin:
            continue
        sig = h_sig.Integral(iStart, iEnd)
        bkg = h_bkg.Integral(iStart, iEnd)
        if (sig+bkg) <= 0 or sig_tot <= 0:
            significance = 0
            sig_loss = 0
        else:
            significance = sig / (sig + bkg)**0.5
            sig_loss = sig / sig_tot
        #significance = sig / (sig + bkg)**0.5
        #sig_loss = sig / sig_tot
        print(f"{i} {center:.5f} {significance:.5f} {sig_loss:.5f}")
        x.append(center)
        y.append(significance)
        l.append(sig_loss)


    graph_sign = ROOT.TGraph(len(x), array.array('d', x), array.array('d', y))
    graph_l = ROOT.TGraph(len(x), array.array('d', x), array.array('d', l))
    max_y = max(y)
    max_index = y.index(max_y)
    max_x = x[max_index]
    max_l = l[max_index]

    canvas = ROOT.TCanvas("", "", 800, 800)
    graph_sign.SetMarkerStyle(20)
    graph_sign.SetMarkerColor(ROOT.kBlue)
    graph_sign.GetXaxis().SetRangeUser(xMin, xMax)
    graph_sign.Draw("AP")
    canvas.Update()

    # Add a marker for the maximum point
    max_marker = ROOT.TMarker(max_x, max_y, 20)
    max_marker.SetMarkerColor(ROOT.kRed)
    max_marker.SetMarkerSize(1.5)
    max_marker.Draw()



    rightmax = 1.0
    print(ROOT.gPad.GetUymin())
    print(ROOT.gPad.GetUymax())
    scale =  ROOT.gPad.GetUymax()/rightmax
    rightmin = ROOT.gPad.GetUymin()/ROOT.gPad.GetUymax()
    graph_l.Scale(scale)
    graph_l.SetLineColor(ROOT.kRed)
    graph_l.SetLineWidth(2)
    graph_l.Draw("SAME L")


    axis_r = ROOT.TGaxis(ROOT.gPad.GetUxmax(),ROOT.gPad.GetUymin(), ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymax(),rightmin,rightmax,510,"+L")
    axis_r.SetLineColor(ROOT.kRed)
    axis_r.SetLabelColor(ROOT.kRed)
    axis_r.Draw()


    # Add a text box to indicate the maximum value
    text = ROOT.TLatex()
    text.SetTextSize(0.03)
    text.SetTextColor(ROOT.kBlack)
    text.DrawLatexNDC(0.1, 0.92, f"Max: x = {max_x}, y = {max_y:.5f}, signal loss = {max_l:.5f}")


    suffix = "_reverse" if reverse else ""
    canvas.SaveAs(f"{outDir}/significance/{hName}{suffix}.png")
    canvas.SaveAs(f"{outDir}/significance/{hName}{suffix}.pdf")
    canvas.Close()



if __name__ == "__main__":
    ecm = args.ecm
    cat = args.cat
    lumi = "10.8" if ecm == 240 else "3"

    inputDir = f"output/h_za/histmaker/ecm{ecm}/"
    baseOutDir = f"/home/submit/jaeyserm/public_html/fccee/h_za/plots_{cat}_ecm{ecm}/"
    outDir = baseOutDir


    procs_cfg_240 = {
        "Za"        : ['wzp6_ee_nunuH_HZa_ecm240', 'wzp6_ee_eeH_HZa_ecm240', 'wzp6_ee_tautauH_HZa_ecm240', 'wzp6_ee_ccH_HZa_ecm240', 'wzp6_ee_bbH_HZa_ecm240', 'wzp6_ee_qqH_HZa_ecm240', 'wzp6_ee_ssH_HZa_ecm240', 'wzp6_ee_mumuH_HZa_ecm240'],
        "qqvvZa"    : ['wzp6_ee_nunuH_HZa_ecm240', 'wzp6_ee_ccH_HZa_ecm240', 'wzp6_ee_bbH_HZa_ecm240', 'wzp6_ee_qqH_HZa_ecm240', 'wzp6_ee_ssH_HZa_ecm240'],
        "qqZa"      : ['wzp6_ee_ccH_HZa_ecm240', 'wzp6_ee_bbH_HZa_ecm240', 'wzp6_ee_qqH_HZa_ecm240', 'wzp6_ee_ssH_HZa_ecm240'],
        "vvZa"      : ['wzp6_ee_nunuH_HZa_ecm240'],
        "gaga"      : ['wz3p6_ee_gammagamma_ecm240'],
        "ZZ"      : ['p8_ee_ZZ_ecm240'],
        "Zgamma"    : ['wz3p6_ee_uu_ecm240', 'wz3p6_ee_dd_ecm240', 'wz3p6_ee_cc_ecm240', 'wz3p6_ee_ss_ecm240', 'wz3p6_ee_bb_ecm240', 'wz3p6_ee_tautau_ecm240', 'wz3p6_ee_mumu_ecm240', 'wz3p6_ee_nunu_ecm240', 'wz3p6_ee_ee_Mee_30_150_ecm240'],
        "ZgammaTau"    : ['wz3p6_ee_tautau_ecm240'],
        "ZgammaHadr"    : ['wz3p6_ee_uu_ecm240', 'wz3p6_ee_dd_ecm240', 'wz3p6_ee_cc_ecm240', 'wz3p6_ee_ss_ecm240', 'wz3p6_ee_bb_ecm240'],
        "ZgammaLept"    : ['wz3p6_ee_mumu_ecm240', 'wz3p6_ee_nunu_ecm240', 'wz3p6_ee_ee_Mee_30_150_ecm240'],

    }

    procs_cfg_365 = {
        #"Za"        : ['wzp6_ee_nunuH_HZa_ecm365', 'wzp6_ee_eeH_HZa_ecm365', 'wzp6_ee_tautauH_HZa_ecm365', 'wzp6_ee_ccH_HZa_ecm365', 'wzp6_ee_bbH_HZa_ecm365', 'wzp6_ee_qqH_HZa_ecm365', 'wzp6_ee_ssH_HZa_ecm365', 'wzp6_ee_mumuH_HZa_ecm365'],
        "Za"        : ['wzp6_ee_nunuH_HZa_ecm365', 'wzp6_ee_ccH_HZa_ecm365', 'wzp6_ee_bbH_HZa_ecm365', 'wzp6_ee_qqH_HZa_ecm365', 'wzp6_ee_ssH_HZa_ecm365'],
        "qqZa"        : ['wzp6_ee_ccH_HZa_ecm365', 'wzp6_ee_bbH_HZa_ecm365', 'wzp6_ee_qqH_HZa_ecm365', 'wzp6_ee_ssH_HZa_ecm365'],
        "vvZa"        : ['wzp6_ee_nunuH_HZa_ecm365'],
        "ZHZa"        : ['wzp6_ee_numunumuH_HZa_ecm365', 'wzp6_ee_ccH_HZa_ecm365', 'wzp6_ee_bbH_HZa_ecm365', 'wzp6_ee_qqH_HZa_ecm365', 'wzp6_ee_ssH_HZa_ecm365'],
        "VBFZa"         : ['wzp6_ee_nuenueVBFH_HZa_ecm365'],
        "gaga"      : ['wz3p6_ee_gammagamma_ecm365'],
        "ZZ"      : ['p8_ee_ZZ_ecm365'],
        "Zgamma"    : ['wz3p6_ee_uu_ecm365', 'wz3p6_ee_dd_ecm365', 'wz3p6_ee_cc_ecm365', 'wz3p6_ee_ss_ecm365', 'wz3p6_ee_bb_ecm365', 'wz3p6_ee_tautau_ecm365', 'wz3p6_ee_mumu_ecm365', 'wz3p6_ee_nunu_ecm365', 'wz3p6_ee_ee_Mee_30_150_ecm365'],
        "ZgammaTau"    : ['wz3p6_ee_tautau_ecm365'],
        "ZgammaHadr"    : ['wz3p6_ee_uu_ecm365', 'wz3p6_ee_dd_ecm365', 'wz3p6_ee_cc_ecm365', 'wz3p6_ee_ss_ecm365', 'wz3p6_ee_bb_ecm365'],
        "ZgammaLept"    : ['wz3p6_ee_mumu_ecm365', 'wz3p6_ee_nunu_ecm365', 'wz3p6_ee_ee_Mee_30_150_ecm365'],

    }

    # colors from https://github.com/mpetroff/accessible-color-cycles
    procs_colors = {
        "Za"            : ROOT.TColor.GetColor("#e42536"),
        "qqvvZa"        : ROOT.TColor.GetColor("#e42536"),
        "vvZa"          : ROOT.TColor.GetColor("#e42536"),
        "ZHZa"        : ROOT.TColor.GetColor("#e42536"),
        "VBFZa"        : ROOT.TColor.GetColor("#e42536"),
        "qqZa"          : ROOT.TColor.GetColor("#e42536"),
        "Zgamma"        : ROOT.TColor.GetColor("#5790fc"),
        "gaga"          : ROOT.TColor.GetColor("#9c9ca1"), #9c9ca1

        "ZgammaTau"     : ROOT.TColor.GetColor("#5790fc"),
        "ZgammaHadr"    : ROOT.TColor.GetColor("#5790fc"),
        "ZgammaLept"    : ROOT.TColor.GetColor("#964a8b"),
        
        "ZZ"      : ROOT.TColor.GetColor("#f89c20")
    }


    procs_labels = {
        "Za"        : "ZH(Z#gamma)",
        "qqvvZa"    : "Z(#nu#nu/qq)H(Z#gamma)",
        "vvZa"      : "Z(#nu#nu)H(Z#gamma)",
        "ZHZa"      : "Z(#nu#nu/qq)H(Z#gamma)",
        "VBFZa"      : "#nu#nuH(Z#gamma)",
        "qqZa"      : "Z(qq)H(Z#gamma)",
        "Zgamma"    : "Z/#gamma^{*} #rightarrow f#bar{f}+#gamma(#gamma)",
        "gaga"      : "#gamma#gamma",
        
        "ZgammaTau"     : "Z/#gamma^{*} #rightarrow taus",
        "ZgammaHadr"    : "Z/#gamma^{*} #rightarrow hadrons",
        "ZgammaLept"    : "Z/#gamma^{*} #rightarrow leptons",
        "ZZ"            : "ZZ",
    }

    if ecm == 240:
        procs_cfg = procs_cfg_240
    else:
        procs_cfg = procs_cfg_365

    
    if cat == "qqvv":

        #procs = ["Za", "Zgamma", "gaga"]
        #procs = ["Za", "ZZ", "ZgammaHadr", "ZgammaTau", "ZgammaLept", "gaga"]
        procs = ["qqvvZa", "ZZ", "Zgamma", "gaga"] # 240
        procs = ["ZHZa", "ZZ", "Zgamma", "gaga"] # 365
        procs = ["VBFZa", "ZZ", "Zgamma", "gaga"] # 365

        # cutflow
        cuts = ["cut0", "cut1", "cut2", "cut3", "cut4", "cut5", "cut6", "cut7", "cut8", "cut9", "cut10", "cut11"]
        labels = ["All events", "#geq 1 #gamma", "1 #gamma", "#gamma iso", "#gamma cos(#theta) < 0.85", "cos#theta_{miss} < 0.995", "Lepton veto", "Missing mass", "d_{12}", "2 jets", "Jet ID", "m_{qq}"]
        #makeCutFlow()
        makeCutFlow(hName="qqvv_cutFlow", cuts=cuts, labels=labels, sig_scale=100, yMin=1e1, yMax=1e10)

        # significance
        if True:

            #significance("photons_veto_n", 0, 10, reverse=True)
            #significance("photons_veto_costheta", 0, 1, reverse=True)
            
            significance("leading_photon_p_nOne", 0, 100)
            significance("leading_photon_p_nOne", 0, 100, reverse=True)
            
            
            significance("photon_iso_nOne", 0, 1.5)
            significance("photon_iso_nOne", 0, 1.5, reverse=True)
            significance("photon_costheta_nOne", 0.5, 1, reverse=True)
            significance("cosThetaEmiss_nOne", 0.9, 1, reverse=True)
            
            # vv missing mass optimization
            significance("missingMass_nOne", 50, 150)
            significance("missingMass_nOne", 50, 150, reverse=True)
            
            significance("qqvv_dmerge_12", 0, 5000)
            
            significance("qqvv_qq_m_nOne", 50, 150)
            significance("qqvv_qq_m_nOne", 50, 150, reverse=True)
            
            #significance("qqvv_mva_score_split_chi2", 0, 1)
            #significance("qqvv_mva_score_split_chi2", 0, 1, reverse=True)




        outDir = f"{baseOutDir}/"
        makePlot("leading_photon_p_nOne", xMin=0, xMax=100, yMin=1e-1, yMax=-1, xLabel="Leading photon momentum (GeV)", logY=True)
        makePlot("photons_n_nOne", xMin=0, xMax=5, yMin=1e-1, yMax=-1, xLabel="Photon multiplicity", logY=True)
        makePlot("photon_costheta_nOne", xMin=0, xMax=1, yMin=1e-1, yMax=-1, xLabel="cos#theta photon", logY=True, rebin=10)
        makePlot("photon_iso_nOne", xMin=0, xMax=5, yMin=1e-1, yMax=-1, xLabel="Relative photon isolation", logY=True)
        makePlot("cosThetaEmiss_nOne", xMin=0.9, xMax=1, yMin=1e-1, yMax=-1, xLabel="cos#theta_{miss}", logY=True)
        makePlot("missingMass_nOne", xMin=0, xMax=220 if ecm==240 else 400, yMin=1e-1, yMax=-1, xLabel="Missing mass (GeV)",logY=True)


        makePlot("qqvv_dmerge_01", xMin=0, xMax=1000, yMin=1e-1, yMax=-1, xLabel="d_{01}", logY=True, rebin=100)
        makePlot("qqvv_dmerge_12", xMin=0, xMax=15000, yMin=1e-1, yMax=-1, xLabel="d_{12}", logY=True, rebin=100)
        makePlot("qqvv_dmerge_23", xMin=0, xMax=5000, yMin=1e-1, yMax=-1, xLabel="d_{23}", logY=True, rebin=100)


        makePlot("qqvv_jets_n_nOne", xMin=0, xMax=5, yMin=1e-1, yMax=-1, xLabel="jets_n_nOne", logY=True)
        makePlot("qqvv_jet1_p_nOne", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="Leading jet momentum (GeV)", logY=True)
        makePlot("qqvv_jet2_p_nOne", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="Subleading jet momentum (GeV)", logY=True)
        makePlot("qqvv_jet1_nc_nOne", xMin=0, xMax=75, yMin=1e-1, yMax=-1, xLabel="Leading jet constituents", logY=True)
        makePlot("qqvv_jet2_nc_nOne", xMin=0, xMax=75, yMin=1e-1, yMax=-1, xLabel="Subleading jet constituents", logY=True)


        makePlot("qqvv_qq_m_nOne", xMin=40, xMax=150, yMin=1e-1, yMax=-1, xLabel="m_{qq} (GeV)", logY=True)
        makePlot("qqvv_qq_costheta_nOne", xMin=0, xMax=1, yMin=1e-1, yMax=-1, xLabel="qq_costheta_nOne", logY=True)
        makePlot("qqvv_qqa_recoil_m_nOne", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="qqa_recoil_m_nOne", logY=True)
        makePlot("qqvv_vva_recoil_m_nOne", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="qqvv_vva_recoil_m_nOne", logY=True)


        makePlot("qqvv_qq_p", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="qq_p", logY=True)
        makePlot("qqvv_qqa_m", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="qqa invariant mass (GeV)", logY=True)
        makePlot("qqvv_qqa_p", xMin=0, xMax=100, yMin=1e-1, yMax=-1, xLabel="qqa momentum (GeV)", logY=True)
        makePlot("qqvv_qq_recoil_m", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="qq_recoil_m", logY=True)
        
        makePlot("qqvv_vva_m", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="vva_m", logY=True)
        makePlot("qqvv_vva_p", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="vva_p", logY=True)
        makePlot("qqvv_vv_recoil_m",  xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="vv_recoil_m", logY=True)

        makePlot("qqvv_vv_trans",  xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="qqvv_vv_trans", logY=True)
        makePlot("qqvv_qq_trans",  xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="qqvv_qq_trans", logY=True)
        makePlot("qqvv_vv_long",  xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="qqvv_vv_long", logY=True)
        makePlot("qqvv_qq_long",  xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="qqvv_qq_long", logY=True)


        # all photon plots
        outDir = f"{baseOutDir}/photons_all"
        makePlot("photons_all_p", xMin=0, xMax=100, yMin=1e-1, yMax=-1, xLabel="photons_all_p", logY=True)
        makePlot("photons_all_n", xMin=0, xMax=10, yMin=1e-1, yMax=-1, xLabel="photons_all_n", logY=True)
        makePlot("photons_all_costheta", xMin=0, xMax=1, yMin=1e-1, yMax=-1, xLabel="photons_all_costheta", logY=True)

        # photon plots
        outDir = f"{baseOutDir}/photons"
        makePlot("photons_p", xMin=0, xMax=100, yMin=1e-1, yMax=-1, xLabel="photons_p", logY=True)
        makePlot("photons_n", xMin=0, xMax=10, yMin=1e-1, yMax=-1, xLabel="photons_n", logY=True)
        makePlot("photons_costheta", xMin=0, xMax=1, yMin=1e-1, yMax=-1, xLabel="photons_costheta", logY=True)

        # veto photon plots
        outDir = f"{baseOutDir}/photons_veto"
        makePlot("photons_veto_p", xMin=0, xMax=100, yMin=1e-1, yMax=-1, xLabel="photons_veto_p", logY=True)
        makePlot("photons_veto_n", xMin=0, xMax=10, yMin=1e-1, yMax=-1, xLabel="photons_veto_n", logY=True)
        makePlot("photons_veto_costheta", xMin=0.8, xMax=1, yMin=1e-3, yMax=-1, xLabel="photons_veto_costheta", logY=True, rebin=10)




        ### chi2 pairing
        outDir = f"{baseOutDir}/pairing_chi2/"
        makePlot("qqvv_H_m", xMin=100, xMax=150, yMin=1e-1, yMax=-1, xLabel="Higgs candidate invariant mass (GeV)", logY=True)
        makePlot("qqvv_H_p", xMin=0, xMax=100, yMin=1e-1, yMax=-1, xLabel="Higgs candidate momentum (GeV)", logY=True)
        makePlot("qqvv_Z_m", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="Associated Z candidate invariant mass (GeV)", logY=True)
        makePlot("qqvv_Z_p", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="Associated Z candidate momentum (GeV)", logY=True)
        makePlot("qqvv_ZH_m", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="qqvv_ZH_m", logY=True)
        makePlot("qqvv_ZH_p", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="qqvv_ZH_p", logY=True)
        makePlot("qqvv_H_recoil_m", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="Higgs candidate recoil mass (GeV)", logY=True)
        makePlot("qqvv_H_recoil_p", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="qqvv_H_recoil_p", logY=True)
        makePlot("qqvv_Z_recoil_m", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="Associated Z candidate recoil mass (GeV)", logY=True)
        makePlot("qqvv_Z_recoil_p", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="qqvv_Z_recoil_p", logY=True)
        makePlot("qqvv_mass_difference_H_Z", xMin=-50, xMax=100, yMin=1e-1, yMax=-1, xLabel="qqvv_mass_difference_H_Z", logY=True)
        makePlot("qqvv_mva_score_split_chi2", xMin=0, xMax=1, yMin=1e-1, yMax=-1, xLabel="MVA score", logY=True, rebin=10)


        quit()
        ### chi2 splitting
        outDir = f"{baseOutDir}/splitting_chi2/hqqa"
        makePlot("qqvv_hqqa_qqa_m_chi2", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="qqvv_hqqa_qqa_m_chi2", logY=True)
        makePlot("qqvv_hqqa_qqa_p_chi2", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="qqvv_hqqa_qqa_p_chi2", logY=True)
        makePlot("qqvv_hqqa_qqa_recoil_m_chi2", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="qqvv_hqqa_qqa_recoil_m_chi2", logY=True)
        makePlot("qqvv_hqqa_final_chi2", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="qqvv_hqqa_final_chi2", logY=True)
        makePlot("qqvv_hqqa_mass_difference_qqa_vv_chi2", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="qqvv_hqqa_mass_difference_qqa_vv_chi2", logY=True)
        makePlot("qqvv_hqqa_mass_difference_vva_qq_chi2", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="qqvv_hqqa_mass_difference_vva_qq_chi2", logY=True)

        outDir = f"{baseOutDir}/splitting_chi2/hvva"
        makePlot("qqvv_hvva_qq_p_chi2", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="qqvv_hvva_qq_p_chi2", logY=True)
        makePlot("qqvv_hvva_qq_recoil_m_chi2", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="qqvv_hvva_qq_recoil_m_chi2", logY=True)
        makePlot("qqvv_hvva_final_chi2", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="qqvv_hvva_final_chi2", logY=True)
        makePlot("qqvv_hvva_mass_difference_qqa_vv_chi2", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="qqvv_hqqa_mass_difference_qqa_vv_chi2", logY=True)
        makePlot("qqvv_hvva_mass_difference_vva_qq_chi2", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="qqvv_hqqa_mass_difference_vva_qq_chi2", logY=True)


        ### MVA splitting
        outDir = f"{baseOutDir}/splitting_mva/hqqa"

        makePlot("qqvv_hqqa_qq_m_mva", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="hqqa_qq_m_mva", logY=True)
        makePlot("qqvv_hqqa_qq_p_mva", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="hqqa_qq_p_mva", logY=True)
        makePlot("qqvv_hqqa_qqa_m_mva", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="hqqa_qqa_m_mva", logY=True)
        makePlot("qqvv_hqqa_qqa_p_mva", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="hqqa_qqa_p_mva", logY=True)
        makePlot("qqvv_hqqa_qq_recoil_m_mva", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="hqqa_qq_recoil_m_mva", logY=True)
        makePlot("qqvv_hqqa_qqa_recoil_za_m_mva", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="hqqa_qqa_recoil_za_m_mva", logY=True)

        makePlot("qqvv_hqqa_vv_m_mva", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="hqqa_vv_m_mva", logY=True)
        makePlot("qqvv_hqqa_vva_m_mva", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="hqqa_vva_m_mva", logY=True)
        makePlot("qqvv_hqqa_vva_p_mva", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="hqqa_vva_p_mva", logY=True)
        makePlot("qqvv_hqqa_vv_recoil_m_mva", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="hqqa_vv_recoil_m_mva", logY=True)

        makePlot("qqvv_hqqa_mass_difference_qqa_vv_mva", xMin=-50, xMax=100, yMin=1e-1, yMax=-1, xLabel="hqqa_mass_difference_qqa_vv",logY=True)
        makePlot("qqvv_hqqa_mass_difference_vva_qq_mva", xMin=-50, xMax=100, yMin=1e-1, yMax=-1, xLabel="hqqa_mass_difference_vva_qq",logY=True)


        makePlot("qqvv_hqqa_dr_vv_qqa_nOne", xMin=0, xMax=5, yMin=1e-1, yMax=-1, xLabel="hqqa_dr_vv_qqa_nOne", logY=True)
        makePlot("qqvv_hqqa_dr_qq_vva_nOne", xMin=0, xMax=5, yMin=1e-1, yMax=-1, xLabel="hqqa_dr_qq_vva_nOne", logY=True)
        makePlot("qqvv_hqqa_dr_a_qq_nOne", xMin=0, xMax=5, yMin=1e-1, yMax=-1, xLabel="hqqa_dr_a_qq_nOne", logY=True)
        makePlot("qqvv_hqqa_dr_a_vv_nOne", xMin=0, xMax=5, yMin=1e-1, yMax=-1, xLabel="hqqa_dr_a_vv_nOne", logY=True)

        makePlot("qqvv_hqqa_mva_score", xMin=0, xMax=1, yMin=1e-2, yMax=-1, xLabel="hqqa_mva_score", logY=True)
        makePlot("qqvv_hqqa_final_1D_mva", xMin=110, xMax=140, yMin=1e-1, yMax=-1, xLabel="hqqa_final_1D_mva", logY=True)





        outDir = f"{baseOutDir}/splitting_mva/hvva"
        makePlot("qqvv_hvva_qq_m_mva", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="hvva_qq_m_mva", logY=True)
        makePlot("qqvv_hvva_qq_p_mva", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="hvva_qq_p_mva", logY=True)
        makePlot("qqvv_hvva_qqa_m_mva", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="hvva_qqa_m_mva", logY=True)
        makePlot("qqvv_hvva_qqa_p_mva", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="hvva_qqa_p_mva", logY=True)
        makePlot("qqvv_hvva_qq_recoil_m_mva", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="hvva_qq_recoil_m_mva", logY=True)
        makePlot("qqvv_hvva_qqa_recoil_za_m_mva", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="hvva_qqa_recoil_za_m_mva", logY=True)

        makePlot("qqvv_hvva_vv_m_mva", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="hvva_vv_m_mva", logY=True)
        makePlot("qqvv_hvva_vva_m_mva", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="hvva_vva_m_mva", logY=True)
        makePlot("qqvv_hvva_vva_p_mva", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="hvva_vva_p_mva", logY=True)
        makePlot("qqvv_hvva_vv_recoil_m_mva", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="hvva_vv_recoil_m_mva", logY=True)

        makePlot("qqvv_hvva_mass_difference_qqa_vv_mva", xMin=-50, xMax=100, yMin=1e-1, yMax=-1, xLabel="hvva_mass_difference_qqa_vv",logY=True)
        makePlot("qqvv_hvva_mass_difference_vva_qq_mva", xMin=-50, xMax=100, yMin=1e-1, yMax=-1, xLabel="",logY=True)

        makePlot("qqvv_hvva_dr_vv_qqa_nOne", xMin=0, xMax=5, yMin=1e-1, yMax=-1, xLabel="hvva_dr_vv_qqa_nOne", logY=True)
        makePlot("qqvv_hvva_dr_qq_vva_nOne", xMin=0, xMax=5, yMin=1e-1, yMax=-1, xLabel="hvva_dr_qq_vva_nOne", logY=True)
        makePlot("qqvv_hvva_dr_a_qq_nOne", xMin=0, xMax=5, yMin=1e-1, yMax=-1, xLabel="hvva_dr_a_qq_nOne", logY=True)
        makePlot("qqvv_hvva_dr_a_vv_nOne", xMin=0, xMax=5, yMin=1e-1, yMax=-1, xLabel="hvva_dr_a_vv_nOne", logY=True)

        makePlot("qqvv_hvva_mva_score", xMin=0, xMax=1, yMin=1e-2, yMax=-1, xLabel="hvva_mva_score", logY=True)
        makePlot("qqvv_hvva_final_1D_mva", xMin=110, xMax=140, yMin=1e-1, yMax=-1, xLabel="hvva_final_1D_mva", logY=True)



    if cat == "vbf":

        inputDir = f"output/h_za/histmaker/ecm{ecm}_vbf/"
        baseOutDir = f"/home/submit/jaeyserm/public_html/fccee/h_za/plots_{cat}_ecm{ecm}/"
        outDir = baseOutDir

        procs = ["VBFZa", "ZZ", "Zgamma", "gaga"] # 365

        # cutflow
        cuts = ["cut0", "cut1", "cut2", "cut3", "cut4", "cut5", "cut6", "cut7", "cut8", "cut9", "cut10", "cut11"]
        labels = ["All events", "#geq 1 #gamma", "#gamma iso", "#gamma cos(#theta) < 0.85", "cos#theta_{miss} < 0.995", "Lepton veto", "Missing mass", "d_{12}", "2 jets", "Jet ID", "m_{qq}"]
        makeCutFlow(hName="cutFlow", cuts=cuts, labels=labels, sig_scale=100, yMin=1e1, yMax=1e10)

        # significance
        if True:

            #significance("photons_veto_n", 0, 10, reverse=True)
            #significance("photons_veto_costheta", 0, 1, reverse=True)
            
            #significance("leading_photon_p_nOne", 0, 100)
            #significance("leading_photon_p_nOne", 0, 100, reverse=True)
            
            
            #significance("photon_iso_nOne", 0, 1.5)
            significance("photon_iso_nOne", 0, 1.5, reverse=True)
            significance("photon_costheta_nOne", 0.5, 1, reverse=True)
            significance("cosThetaEmiss_nOne", 0.9, 1, reverse=True)
            
            significance("dmerge_12", 0, 5000)
            
            significance("cos_qqa_nOne", 0.5, 1, reverse=True)
            
            significance("vv_trans_nOne", 0, 50)
            
            # vv missing mass optimization
            #significance("missingMass_nOne", 50, 150)
            #significance("missingMass_nOne", 50, 150, reverse=True)
            
            #
            
            #significance("qqvv_qq_m_nOne", 50, 150)
            #significance("qqvv_qq_m_nOne", 50, 150, reverse=True)
            
            #significance("qqvv_mva_score_split_chi2", 0, 1)
            #significance("qqvv_mva_score_split_chi2", 0, 1, reverse=True)




        outDir = f"{baseOutDir}/"
        makePlot("leading_photon_p_nOne", xMin=0, xMax=100, yMin=1e-2, yMax=-1, xLabel="Leading photon momentum (GeV)", logY=True)
        makePlot("photons_n_nOne", xMin=0, xMax=5, yMin=1e-2, yMax=-1, xLabel="Photon multiplicity", logY=True)
        makePlot("photon_costheta_nOne", xMin=0, xMax=1, yMin=1e-2, yMax=-1, xLabel="cos#theta photon", logY=True, rebin=10)
        makePlot("photon_iso_nOne", xMin=0, xMax=5, yMin=1e-2, yMax=-1, xLabel="Relative photon isolation", logY=True)
        makePlot("cosThetaEmiss_nOne", xMin=0.9, xMax=1, yMin=1e-2, yMax=-1, xLabel="cos#theta_{miss}", logY=True)
        makePlot("missingMass_nOne", xMin=0, xMax=220 if ecm==240 else 400, yMin=1e-1, yMax=-1, xLabel="Missing mass (GeV)",logY=True)


        makePlot("dmerge_01", xMin=0, xMax=1000, yMin=1e-2, yMax=-1, xLabel="d_{01}", logY=True, rebin=100)
        makePlot("dmerge_12", xMin=0, xMax=15000, yMin=1e-2, yMax=-1, xLabel="d_{12}", logY=True, rebin=100)
        makePlot("dmerge_23", xMin=0, xMax=5000, yMin=1e2, yMax=-1, xLabel="d_{23}", logY=True, rebin=100)


        makePlot("jets_n_nOne", xMin=0, xMax=5, yMin=1e-2, yMax=-1, xLabel="jets_n_nOne", logY=True)
        makePlot("jet1_p_nOne", xMin=0, xMax=250, yMin=1e-2, yMax=-1, xLabel="Leading jet momentum (GeV)", logY=True)
        makePlot("jet2_p_nOne", xMin=0, xMax=250, yMin=1e-2, yMax=-1, xLabel="Subleading jet momentum (GeV)", logY=True)
        makePlot("jet1_nc_nOne", xMin=0, xMax=75, yMin=1e-2, yMax=-1, xLabel="Leading jet constituents", logY=True)
        makePlot("jet2_nc_nOne", xMin=0, xMax=75, yMin=1e-2, yMax=-1, xLabel="Subleading jet constituents", logY=True)


        #makePlot("qq_m", xMin=40, xMax=150, yMin=1e-2, yMax=-1, xLabel="m_{qq} (GeV)", logY=True)
        #makePlot("qq_costheta_nOne", xMin=0, xMax=1, yMin=1e-2, yMax=-1, xLabel="qq_costheta_nOne", logY=True)


        makePlot("qq_p", xMin=0, xMax=250, yMin=1e-2, yMax=-1, xLabel="qq_p", logY=True)
        makePlot("qqa_m", xMin=0, xMax=250, yMin=1e-2, yMax=-1, xLabel="qqa invariant mass (GeV)", logY=True)
        makePlot("qqa_p", xMin=0, xMax=100, yMin=1e-2, yMax=-1, xLabel="qqa momentum (GeV)", logY=True)


        makePlot("vv_trans_nOne", xMin=0, xMax=250, yMin=1e-2, yMax=-1, xLabel="vv_trans", logY=True)
        makePlot("qq_trans", xMin=0, xMax=250, yMin=1e-2, yMax=-1, xLabel="qq_trans", logY=True)
        makePlot("vv_long", xMin=0, xMax=100, yMin=1e-2, yMax=-1, xLabel="vv_long", logY=True)
        makePlot("qq_long", xMin=0, xMax=100, yMin=1e-2, yMax=-1, xLabel="qq_long", logY=True)


        makePlot("dr_vv_qqa", xMin=0, xMax=5, yMin=1e-1, yMax=-1, xLabel="dr_vv_qqa", logY=True)
        makePlot("dr_vv_qq", xMin=0, xMax=5, yMin=1e-1, yMax=-1, xLabel="dr_vv_qq", logY=True)
        makePlot("dr_a_qq", xMin=0, xMax=5, yMin=1e-1, yMax=-1, xLabel="dr_a_qq", logY=True)
        makePlot("dr_a_vv", xMin=0, xMax=5, yMin=1e-1, yMax=-1, xLabel="dr_a_vv", logY=True)
       
        makePlot("cos_qqa_nOne", xMin=0, xMax=1, yMin=1e-1, yMax=-1, xLabel="cos_qqa", logY=True, rebin=10)
        makePlot("cos_qq", xMin=0, xMax=1, yMin=1e-1, yMax=-1, xLabel="cos_qq", logY=True, rebin=10)


        quit()

        makePlot("vbf_vv_trans",  xMin=0, xMax=250, yMin=1e-2, yMax=-1, xLabel="vbf_vv_trans", logY=True)
        makePlot("vbf_qq_trans",  xMin=0, xMax=250, yMin=1e-2, yMax=-1, xLabel="vbf_qq_trans", logY=True)
        makePlot("vbf_vv_long",  xMin=0, xMax=250, yMin=1e-2, yMax=-1, xLabel="vbf_vv_long", logY=True)
        makePlot("vbf_qq_long",  xMin=0, xMax=250, yMin=1e-2, yMax=-1, xLabel="vbf_qq_long", logY=True)


        # all photon plots
        outDir = f"{baseOutDir}/photons_all"
        makePlot("photons_all_p", xMin=0, xMax=100, yMin=1e-1, yMax=-1, xLabel="photons_all_p", logY=True)
        makePlot("photons_all_n", xMin=0, xMax=10, yMin=1e-1, yMax=-1, xLabel="photons_all_n", logY=True)
        makePlot("photons_all_costheta", xMin=0, xMax=1, yMin=1e-1, yMax=-1, xLabel="photons_all_costheta", logY=True)

        # photon plots
        outDir = f"{baseOutDir}/photons"
        makePlot("photons_p", xMin=0, xMax=100, yMin=1e-1, yMax=-1, xLabel="photons_p", logY=True)
        makePlot("photons_n", xMin=0, xMax=10, yMin=1e-1, yMax=-1, xLabel="photons_n", logY=True)
        makePlot("photons_costheta", xMin=0, xMax=1, yMin=1e-1, yMax=-1, xLabel="photons_costheta", logY=True)

        # veto photon plots
        outDir = f"{baseOutDir}/photons_veto"
        makePlot("photons_veto_p", xMin=0, xMax=100, yMin=1e-1, yMax=-1, xLabel="photons_veto_p", logY=True)
        makePlot("photons_veto_n", xMin=0, xMax=10, yMin=1e-1, yMax=-1, xLabel="photons_veto_n", logY=True)
        makePlot("photons_veto_costheta", xMin=0.8, xMax=1, yMin=1e-3, yMax=-1, xLabel="photons_veto_costheta", logY=True, rebin=10)











    if cat == "qqqq":

        procs = ["Za", "Zgamma", "gaga"]
        procs = ["Za", "ZZ", "ZgammaHadr", "ZgammaTau", "ZgammaLept", "gaga"]

        # cutflow
        cuts = ["cut0", "cut1", "cut2", "cut3", "cut4", "cut5", "cut6", "cut7", "cut8", "cut9", "cut10", "cut11", "cut12", "cut13", "cut14"]
        labels = ["All events", "#geq 1 #gamma", "1 #gamma", "#gamma iso", "#gamma cos(#theta) < 0.85", "cos#theta_{miss} < 0.995", "Lepton veto", "Missing mass", "d_{34} > 575", "Njets = 4", "Jet ID", "chi2<0", "chi2>0"]
        #makeCutFlow()
        makeCutFlow(hName="qqqq_cutFlow", cuts=cuts, labels=labels, sig_scale=100, yMin=1e1, yMax=1e10)


        # significance
        if True:
            
            #significance("photons_veto_n", 0, 10, reverse=True)
            #significance("photons_veto_costheta", 0, 1, reverse=True)
            
            #significance("leading_photon_p_nOne", 0, 100)
            #significance("leading_photon_p_nOne", 0, 100, reverse=True)
            
            
            #significance("photon_iso_nOne", 0, 1.5)
            #significance("photon_iso_nOne", 0, 1.5, reverse=True)
            #significance("photon_costheta_nOne", 0.5, 1, reverse=True)
            #significance("cosThetaEmiss_nOne", 0.9, 1, reverse=True)
            
            # vv missing mass optimization
            significance("missingMass_nOne", 0, 50)
            significance("missingMass_nOne", 0, 50, reverse=True)
            
            significance("qqqq_jet1_p_nOne", 0, 150)
            significance("qqqq_jet1_p_nOne", 0, 150, reverse=True)
            
            significance("qqqq_jet2_p_nOne", 0, 150)
            significance("qqqq_jet2_p_nOne", 0, 150, reverse=True)
            
            significance("qqqq_jet3_p_nOne", 0, 150)
            significance("qqqq_jet3_p_nOne", 0, 150, reverse=True)
            
            significance("qqqq_jet4_p_nOne", 0, 150)
            significance("qqqq_jet4_p_nOne", 0, 150, reverse=True)
            
            significance("qqqq_acolinearity_2", 0, 1)
            
            significance("qqqq_H_recoil_m", 0, 150)
            significance("qqqq_H_recoil_m", 0, 150, reverse=True)
            
            
            significance("qqqq_thrust_magn_nOne", 0, 1, reverse=True)
            significance("qqqq_thrust_costheta_nOne", 0, 1, reverse=True)
            
            significance("qqqq_dmerge_34_nOne", 0, 2000)
            significance("qqqq_dmerge_45_nOne", 0, 2000)
            
            
            #significance("qq_m_nOne", 50, 150)
            #significance("qq_m_nOne", 50, 150, reverse=True)
            
            #significance("qqa_recoil_m_nOne", 50, 150)
            #significance("qqa_recoil_m_nOne", 50, 150, reverse=True)
            
            #significance("qqa_m", 50, 150)
            #significance("qqa_m", 50, 150, reverse=True)
            
            
            #significance("qq_costheta_nOne", 0, 1, reverse=True)
            
            #significance("hqqa_mva_score", 0, 1)
            #significance("hvva_mva_score", 0, 1)
            
            
        
        outDir = f"{baseOutDir}/"
        # photons
        #makePlot("leading_photon_p_nOne", "leading_photon_p_nOne", xMin=0, xMax=100, yMin=1e-1, yMax=-1, xLabel="Leading photon momentum (GeV)", logY=True)
        #makePlot("photons_n_nOne", "photons_n_nOne", xMin=0, xMax=5, yMin=1e-1, yMax=-1, xLabel="Photon multiplicity", logY=True)
        #makePlot("photon_costheta_nOne", "leading_photon_costheta_nOne_nOne", xMin=0, xMax=1, yMin=1e-1, yMax=-1, xLabel="cos#theta photon", logY=True)
        #makePlot("photon_iso_nOne", "photon_iso_nOne", xMin=0, xMax=5, yMin=1e-1, yMax=-1, xLabel="Photon isolation", logY=True)
        makePlot("cosThetaEmiss_nOne", "cosThetaEmiss_nOne", xMin=0.9, xMax=1, yMin=1e-1, yMax=-1, xLabel="cos#theta_{miss}", logY=True)
        makePlot("missingMass_nOne", "missingMass_nOne", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="Missing mass (GeV)",logY=True)
        #makePlot("missingEnergy_nOne", "missingEnergy_nOne", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="Missing energy (GeV)",logY=True)
        #makePlot("mva_score", "mva_score", xMin=0, xMax=1, yMin=1e-1, yMax=-1, xLabel="mva_score",logY=True)

        makePlot("qqqq_mass_difference_Z1a_Z2", "qqqq_mass_difference_Z1a_Z2", xMin=-50, xMax=100, yMin=1e-1, yMax=-1, xLabel="mass_difference_Z1a_Z2",logY=True)
        makePlot("qqqq_mass_difference_Z2a_Z1", "qqqq_mass_difference_Z2a_Z1", xMin=-50, xMax=100, yMin=1e-1, yMax=-1, xLabel="mass_difference_Z2a_Z1",logY=True)

        outDir = f"{baseOutDir}/jets/"
        makePlot("qqqq_jets_n_nOne",  xMin=0, xMax=5, yMin=1e-1, yMax=-1, xLabel="jets_n", logY=True)
        makePlot("qqqq_jet1_p_nOne", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="jet1_p", logY=True)
        makePlot("qqqq_jet2_p_nOne", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="jet2_p", logY=True)
        makePlot("qqqq_jet3_p_nOne", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="jet3_p", logY=True)
        makePlot("qqqq_jet4_p_nOne", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="jet4_p", logY=True)
        
        makePlot("qqqq_jet1_nc_nOne", xMin=0, xMax=75, yMin=1e-1, yMax=-1, xLabel="jet1_nc", logY=True)
        makePlot("qqqq_jet2_nc_nOne", xMin=0, xMax=75, yMin=1e-1, yMax=-1, xLabel="jet2_nc", logY=True)
        makePlot("qqqq_jet3_nc_nOne", xMin=0, xMax=75, yMin=1e-1, yMax=-1, xLabel="jet3_nc", logY=True)
        makePlot("qqqq_jet4_nc_nOne", xMin=0, xMax=75, yMin=1e-1, yMax=-1, xLabel="jet4_nc", logY=True)


        



        outDir = f"{baseOutDir}/"
        
        makePlot("qqqq_dmerge_01", xMin=0, xMax=50000, yMin=1e-1, yMax=-1, xLabel="dmerge_01", logY=True)
        makePlot("qqqq_dmerge_12", xMin=0, xMax=50000, yMin=1e-1, yMax=-1, xLabel="dmerge_12", logY=True, rebin=100)
        makePlot("qqqq_dmerge_23", xMin=0, xMax=15000, yMin=1e-1, yMax=-1, xLabel="dmerge_23", logY=True, rebin=10)
        makePlot("qqqq_dmerge_34", xMin=0, xMax=5000, yMin=1e-1, yMax=-1, xLabel="dmerge_34", logY=True, rebin=10)
        makePlot("qqqq_dmerge_45", xMin=0, xMax=2000, yMin=1e-1, yMax=-1, xLabel="dmerge_45", logY=True)
        
        makePlot("qqqq_thrust_magn_nOne", xMin=0, xMax=1, yMin=1e-1, yMax=-1, xLabel="qqqq_thrust_magn_nOne", logY=True)
        makePlot("qqqq_thrust_costheta_nOne", xMin=0, xMax=1, yMin=1e-1, yMax=-1, xLabel="qqqq_thrust_costheta_nOne", logY=True, rebin=100)
        
        
        makePlot("qqqq_acolinearity_2", xMin=0, xMax=4, yMin=1e-1, yMax=-1, xLabel="acolinearity", logY=True)
        makePlot("qqqq_acoplanarity_2", xMin=0, xMax=4, yMin=1e-1, yMax=-1, xLabel="acoplanarity", logY=True)
        makePlot("qqqq_dijet_m", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="qqqq_dijet_m", logY=True)
        makePlot("qqqq_dijet_p", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="qqqq_dijet_p", logY=True)
        
        
        makePlot("jets_tot_m_nOne", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="jets_tot_m", logY=True)
        makePlot("jets_tot_p_nOne", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="jets_tot_p", logY=True)

        makePlot("qqqq_Z1_m", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="qqqq_Z1_m", logY=True)
        makePlot("qqqq_Z2_m", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="qqqq_Z2_m", logY=True)
        makePlot("qqqq_Z1_p", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="qqqq_Z1_p", logY=True)
        makePlot("qqqq_Z2_p", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="qqqq_Z2_p", logY=True)


        makePlot("qqqq_H_m", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="qqqq_H_m", logY=True)
        makePlot("qqqq_H_p", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="qqqq_H_p", logY=True)
        makePlot("qqqq_Z_m", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="qqqq_Z_m", logY=True)
        makePlot("qqqq_Z_p", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="qqqq_Z_p", logY=True)

        makePlot("qqqq_ZH_m", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="qqqq_ZH_m", logY=True)
        makePlot("qqqq_ZH_p", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="qqqq_ZH_p", logY=True)

        makePlot("qqqq_H_recoil_m", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="qqqq_H_recoil_m", logY=True)
        makePlot("qqqq_H_recoil_p", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="qqqq_H_recoil_p", logY=True)
        makePlot("qqqq_Z_recoil_m", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="qqqq_Z_recoil_m", logY=True)
        makePlot("qqqq_Z_recoil_p", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="qqqq_Z_recoil_p", logY=True)

        makePlot("qqqq_mass_difference", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="qqqq_mass_difference", logY=True)
        makePlot("qqqq_H_m_final", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="qqqq_H_m_final", logY=True)
        makePlot("qqqq_Z_recoil_m_final", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="qqqq_Z_recoil_m_final", logY=True)



        outDir = f"{baseOutDir}/chi2_hz1/"
        makePlot("qqqq_H_recoil_m_hz1", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="qqqq_H_recoil_m_hz1", logY=True)
        makePlot("qqqq_H_recoil_p_hz1", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="qqqq_H_recoil_p_hz1", logY=True)
        makePlot("qqqq_Z_recoil_m_hz1", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="qqqq_Z_recoil_m_hz1", logY=True)
        makePlot("qqqq_Z_recoil_p_hz1", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="qqqq_Z_recoil_p_hz1", logY=True)
        
        makePlot("qqqq_H_m_hz1", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="qqqq_H_m_hz1", logY=True)
        makePlot("qqqq_H_p_hz1", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="qqqq_H_p_hz1", logY=True)
        makePlot("qqqq_Z_m_hz1", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="qqqq_Z_m_hz1", logY=True)
        makePlot("qqqq_Z_p_hz1", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="qqqq_Z_p_hz1", logY=True)
        makePlot("qqqq_ZH_m_hz1", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="qqqq_ZH_m_hz1", logY=True)
        makePlot("qqqq_ZH_p_hz1", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="qqqq_ZH_p_hz1", logY=True)


        outDir = f"{baseOutDir}/chi2_hz2/"
        makePlot("qqqq_H_recoil_m_hz2", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="qqqq_H_recoil_m_hz2", logY=True)
        makePlot("qqqq_H_recoil_p_hz2", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="qqqq_H_recoil_p_hz2", logY=True)
        makePlot("qqqq_Z_recoil_m_hz2", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="qqqq_Z_recoil_m_hz2", logY=True)
        makePlot("qqqq_Z_recoil_p_hz2", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="qqqq_Z_recoil_p_hz2", logY=True)
        
        makePlot("qqqq_H_m_hz2", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="qqqq_H_m_hz2", logY=True)
        makePlot("qqqq_H_p_hz2", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="qqqq_H_p_hz2", logY=True)
        makePlot("qqqq_Z_m_hz2", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="qqqq_Z_m_hz2", logY=True)
        makePlot("qqqq_Z_p_hz2", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="qqqq_Z_p_hz2", logY=True)
        makePlot("qqqq_ZH_m_hz2", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="qqqq_ZH_m_hz2", logY=True)
        makePlot("qqqq_ZH_p_hz2", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="qqqq_ZH_p_hz2", logY=True)


        quit()

        makePlot("qq_m_nOne", "qq_m_nOne", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="qq_m_nOne", logY=True)
        makePlot("qq_costheta_nOne", "qq_costheta_nOne", xMin=0, xMax=1, yMin=1e-1, yMax=-1, xLabel="qq_costheta_nOne", logY=True)
        makePlot("qqa_recoil_m_nOne", "qqa_recoil_m_nOne", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="qqa_recoil_m_nOne", logY=True)

        
        makePlot("qq_p", "qq_p_nOne", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="qq_p", logY=True)
        makePlot("qqa_m", "qqa_m_nOne", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="qqa_m", logY=True)
        makePlot("qqa_p", "qqa_p_nOne", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="qqa_p", logY=True)
        makePlot("qq_recoil_m", "qq_recoil_m_nOne", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="qq_recoil_m", logY=True)
        
        makePlot("vva_m", "vva_m_nOne", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="vva_m", logY=True)
        makePlot("vva_p", "vva_p_nOne", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="vva_p", logY=True)
        makePlot("vv_recoil_m", "vv_recoil_m_nOne", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="vv_recoil_m", logY=True)


        #makePlot("hqqa_final_1D_mva", "hqqa_final_1D_mva", xMin=110, xMax=140, yMin=1e-1, yMax=-1, xLabel="hqqa_final_1D_mva", logY=True)
        #makePlot("hvva_final_1D_mva", "hvva_final_1D_mva", xMin=110, xMax=140, yMin=1e-1, yMax=-1, xLabel="hvva_final_1D_mva", logY=True)




        outDir = f"{baseOutDir}/plots/hqqa/"
        makePlot("hqqa_qq_m_mva", "hqqa_qq_m_mva", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="hqqa_qq_m_mva", logY=True)
        makePlot("hqqa_qq_p_mva", "hqqa_qq_p_mva", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="hqqa_qq_p_mva", logY=True)
        makePlot("hqqa_qqa_m_mva", "hqqa_qqa_m_mva", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="hqqa_qqa_m_mva", logY=True)
        makePlot("hqqa_qqa_p_mva", "hqqa_qqa_p_mva", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="hqqa_qqa_p_mva", logY=True)
        makePlot("hqqa_qq_recoil_m_mva", "hqqa_qq_recoil_m_mva", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="hqqa_qq_recoil_m_mva", logY=True)
        makePlot("hqqa_qqa_recoil_za_m_mva", "hqqa_qqa_recoil_za_m_mva", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="hqqa_qqa_recoil_za_m_mva", logY=True)

        makePlot("hqqa_vv_m_mva", "hqqa_vv_m_mva", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="hqqa_vv_m_mva", logY=True)
        makePlot("hqqa_vva_m_mva", "hqqa_vva_m_mva", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="hqqa_vva_m_mva", logY=True)
        makePlot("hqqa_vva_p_mva", "hqqa_vva_p_mva", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="hqqa_vva_p_mva", logY=True)
        makePlot("hqqa_vv_recoil_m_mva", "hqqa_vv_recoil_m_mva", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="hqqa_vv_recoil_m_mva", logY=True)

        makePlot("hqqa_mass_difference_qqa_vv", "hqqa_mass_difference_qqa_vv", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="hqqa_mass_difference_qqa_vv", logY=True)
        makePlot("hqqa_mass_difference_vva_qq", "hqqa_mass_difference_vva_qq", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="hqqa_mass_difference_vva_qq", logY=True)


        makePlot("hqqa_dr_vv_qqa_nOne", "hqqa_dr_vv_qqa_nOne", xMin=0, xMax=5, yMin=1e-1, yMax=-1, xLabel="hqqa_dr_vv_qqa_nOne", logY=True)
        makePlot("hqqa_dr_qq_vva_nOne", "hqqa_dr_qq_vva_nOne", xMin=0, xMax=5, yMin=1e-1, yMax=-1, xLabel="hqqa_dr_qq_vva_nOne", logY=True)
        makePlot("hqqa_dr_a_qq_nOne", "hqqa_dr_a_qq_nOne", xMin=0, xMax=5, yMin=1e-1, yMax=-1, xLabel="hqqa_dr_a_qq_nOne", logY=True)
        makePlot("hqqa_dr_a_vv_nOne", "hqqa_dr_a_vv_nOne", xMin=0, xMax=5, yMin=1e-1, yMax=-1, xLabel="hqqa_dr_a_vv_nOne", logY=True)

        makePlot("hqqa_mva_score", "hqqa_mva_score", xMin=0, xMax=1, yMin=1e-2, yMax=-1, xLabel="hqqa_mva_score", logY=True)



        outDir = f"{baseOutDir}/plots/hvva/"
        makePlot("hvva_qq_m_mva", "hvva_qq_m_mva", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="hvva_qq_m_mva", logY=True)
        makePlot("hvva_qq_p_mva", "hvva_qq_p_mva", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="hvva_qq_p_mva", logY=True)
        makePlot("hvva_qqa_m_mva", "hvva_qqa_m_mva", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="hvva_qqa_m_mva", logY=True)
        makePlot("hvva_qqa_p_mva", "hvva_qqa_p_mva", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="hvva_qqa_p_mva", logY=True)
        makePlot("hvva_qq_recoil_m_mva", "hvva_qq_recoil_m_mva", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="hvva_qq_recoil_m_mva", logY=True)
        makePlot("hvva_qqa_recoil_za_m_mva", "hvva_qqa_recoil_za_m_mva", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="hvva_qqa_recoil_za_m_mva", logY=True)

        makePlot("hvva_vv_m_mva", "hvva_vv_m_mva", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="hvva_vv_m_mva", logY=True)
        makePlot("hvva_vva_m_mva", "hvva_vva_m_mva", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="hvva_vva_m_mva", logY=True)
        makePlot("hvva_vva_p_mva", "hvva_vva_p_mva", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="hvva_vva_p_mva", logY=True)
        makePlot("hvva_vv_recoil_m_mva", "hvva_vv_recoil_m_mva", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="hvva_vv_recoil_m_mva", logY=True)

        makePlot("hvva_mass_difference_qqa_vv", "hvva_mass_difference_qqa_vv", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="hvva_mass_difference_qqa_vv", logY=True)
        makePlot("hvva_mass_difference_vva_qq", "hvva_mass_difference_vva_qq", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="hvva_mass_difference_vva_qq", logY=True)

        makePlot("hvva_dr_vv_qqa_nOne", "hvva_dr_vv_qqa_nOne", xMin=0, xMax=5, yMin=1e-1, yMax=-1, xLabel="hvva_dr_vv_qqa_nOne", logY=True)
        makePlot("hvva_dr_qq_vva_nOne", "hvva_dr_qq_vva_nOne", xMin=0, xMax=5, yMin=1e-1, yMax=-1, xLabel="hvva_dr_qq_vva_nOne", logY=True)
        makePlot("hvva_dr_a_qq_nOne", "hvva_dr_a_qq_nOne", xMin=0, xMax=5, yMin=1e-1, yMax=-1, xLabel="hvva_dr_a_qq_nOne", logY=True)
        makePlot("hvva_dr_a_vv_nOne", "hvva_dr_a_vv_nOne", xMin=0, xMax=5, yMin=1e-1, yMax=-1, xLabel="hvva_dr_a_vv_nOne", logY=True)

        makePlot("hvva_mva_score", "hvva_mva_score", xMin=0, xMax=1, yMin=1e-2, yMax=-1, xLabel="hvva_mva_score", logY=True)