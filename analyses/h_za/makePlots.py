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

def makeCutFlow():

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



def makePlot(hName, outName, xMin=0, xMax=100, yMin=1, yMax=1e5, xLabel="xlabel", yLabel="Events", logX=False, logY=True, rebin=1, legPos=[0.3, 0.75, 0.9, 0.9]):


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


def significance(hName, xMin=-10000, xMax=10000, reverse=False):

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

if __name__ == "__main__":

    inputDir = "output/h_za/histmaker/"
    outDir = "/home/submit/jaeyserm/public_html/fccee/h_za/plots/"
    
    ecm, lumi = 240, "10.8"


    sigs = ['wzp6_ee_nunuH_HZa_ecm240', 'wzp6_ee_eeH_HZa_ecm240', 'wzp6_ee_tautauH_HZa_ecm240', 'wzp6_ee_ccH_HZa_ecm240', 'wzp6_ee_bbH_HZa_ecm240', 'wzp6_ee_qqH_HZa_ecm240', 'wzp6_ee_ssH_HZa_ecm240', 'wzp6_ee_mumuH_HZa_ecm240']
    sigs = ['wzp6_ee_ccH_HZa_ecm240', 'wzp6_ee_bbH_HZa_ecm240', 'wzp6_ee_qqH_HZa_ecm240', 'wzp6_ee_ssH_HZa_ecm240']
    #sigs = ['wzp6_ee_nunuH_HZa_ecm240']
    sig_scale, sig_legend = 100, "ZH(Z#gamma) (#times 100)"

    bkgs = ["gaga", "Zgamma"]

    bkgs_legends = ["#gamma#gamma", "Z/#gamma^{*} #rightarrow f#bar{f}+#gamma(#gamma)"]
    bkgs_colors = [ROOT.TColor.GetColor(248, 206, 104), ROOT.TColor.GetColor(222, 90, 106), ROOT.TColor.GetColor(100, 192, 232), ROOT.TColor.GetColor(155, 152, 204)]
    bgks_cfg = { 
        "gaga"      : ['wzp6_ee_gammagamma_ecm240'],
        "Zgamma"    : ['kkmcee_ee_uu_ecm240', 'kkmcee_ee_dd_ecm240', 'kkmcee_ee_cc_ecm240', 'kkmcee_ee_ss_ecm240', 'kkmcee_ee_bb_ecm240', 'kkmcee_ee_tautau_ecm240', 'kkmcee_ee_mumu_ecm240', 'kkmcee_ee_nuenue_ecm240', 'kkmcee_ee_numunumu_ecm240', 'kkmcee_ee_nutaunutau_ecm240'],
        #"Zgamma"    : ['kkmcee_ee_uu_ecm240', 'kkmcee_ee_dd_ecm240', 'kkmcee_ee_cc_ecm240', 'kkmcee_ee_ss_ecm240', 'kkmcee_ee_bb_ecm240', 'kkmcee_ee_tautau_ecm240', 'kkmcee_ee_mumu_ecm240', 'kkmcee_ee_nuenue_ecm240', 'kkmcee_ee_numunumu_ecm240', 'kkmcee_ee_nutaunutau_ecm240'],
        "ZgammaKKMC"    : ['kkmcee_ee_uu_ecm240', 'kkmcee_ee_dd_ecm240', 'kkmcee_ee_cc_ecm240', 'kkmcee_ee_ss_ecm240', 'kkmcee_ee_bb_ecm240', 'kkmcee_ee_nuenue_ecm240', 'kkmcee_ee_numunumu_ecm240', 'kkmcee_ee_nutaunutau_ecm240', 'kkmcee_ee_tautau_ecm240'],
        "ZgammaP8"  : ['p8_ee_Zqq_ecm240'],
    }

    # cutflow
    cuts = ["cut0", "cut1", "cut2", "cut3", "cut4", "cut5", "cut6", "cut7", "cut8"]
    labels = ["All events", "#geq 1 #gamma", "#gamma p", "#gamma iso", "#gamma costheta", "20 < p_{#gamma#gamma} < 65", "acolinearity > 0.10", "cos#theta_{miss} < 0.99", "Missing energy", "110 < m_{#gamma#gamma} < 130"]
    makeCutFlow()


    sig_scale, sig_legend = 1, "ZH(Z#gamma)"

    # significance
    if True:
        significance("photon_iso", 0, 1, reverse=True)
        significance("photon_costheta", 0.5, 1, reverse=True)
        significance("cosThetaEmiss", 0.9, 1, reverse=True)
        
        # vv missing mass optimization
        significance("missingMass", 50, 150)
        significance("missingMass", 50, 150, reverse=True)
        
        
        
        # vvqq optimization
        significance("vvqq_qq_m", 50, 130)
        significance("vvqq_qq_m", 50, 130, reverse=True)
        significance("vvqq_qq_costheta", 0, 1, reverse=True)
        
        #quit()
        #significance("leading_photon_costheta_nOne", 0.5, 1, reverse=True)
        
        
        
        
        # vv_qq_recoil_za_m
        #significance("vv_qq_recoil_za_m_nOne", 50, 130)
        #significance("vv_qq_recoil_za_m_nOne", 50, 130, reverse=True)
        
        # vv_mumu_recoil_za_m
        #significance("vv_mumu_recoil_za_m", 50, 130)
        #significance("vv_mumu_recoil_za_m", 50, 130, reverse=True)

    # photons
    makePlot("leading_photon_p", "leading_photon_p", xMin=0, xMax=100, yMin=1e-1, yMax=-1, xLabel="Leading photon momentum (GeV)", logY=True)
    makePlot("photons_n", "photons_n", xMin=0, xMax=5, yMin=1e-1, yMax=-1, xLabel="Photon multiplicity", logY=True)
    makePlot("photon_costheta", "leading_photon_costheta_nOne", xMin=0, xMax=1, yMin=1e-1, yMax=-1, xLabel="cos#theta photon", logY=True)
    makePlot("photon_iso", "photon_iso", xMin=0, xMax=5, yMin=1e-1, yMax=-1, xLabel="Photon isolation", logY=True)
    makePlot("cosThetaEmiss", "cosThetaEmiss", xMin=0.9, xMax=1, yMin=1e-1, yMax=-1, xLabel="cos#theta_{miss}", logY=True)
    makePlot("missingEnergy", "missingEnergy", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="Missing energy (GeV)", logY=True)
    makePlot("missingMass", "missingMass", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="Missing mass (GeV)",logY=True)
    
    
    
    
    
    
    
    
    
    
    
    ## comparison plots kkmc-p8
    '''
    outDir_comparison = "/home/submit/jaeyserm/public_html/fccee/h_za/plots_kkmc_p8/"
    comparison("photons_n_nOne", "photons_n_nOne", xMin=0, xMax=10, yMax=-1, xLabel="Photon multiplicity", yLabel="Events", logY=False)
    comparison("leading_photon_p_nOne", "leading_photon_p_nOne", xMin=0, xMax=150, yMax=-1, xLabel="Leading photon momentum (GeV)", yLabel="Events", logY=False)
    comparison("leading_photon_iso_nOne", "leading_photon_iso_nOne", xMin=0, xMax=5, yMax=-1, xLabel="Leading photon isolation", yLabel="Events", logY=False)
    comparison("leading_photon_costheta_nOne", "leading_photon_costheta_nOne", xMin=0, yMin=1e-1, yMax=-1, xLabel="cos#theta leading photon", yLabel="Events", logY=False)
    comparison("cosThetaMiss_nOne", "cosThetaMiss_nOne", xMin=0.9, xMax=1, yMax=-1, xLabel="cos#theta_{miss}", yLabel="Events", logY=False, rebin=1)
    comparison("missingEnergy", "missingEnergy", xMin=0, xMax=250, yMax=-1, xLabel="Missing energy (GeV)", yLabel="Events", logY=False, rebin=1)
    comparison("missingMass", "missingMass", xMin=0, xMax=250, yMax=-1, xLabel="Missing mass (GeV)", yLabel="Events", logY=False, rebin=1)
    '''



    # vv_qq
    outDir = "/home/submit/jaeyserm/public_html/fccee/h_za/plots/vvqq"
    makePlot("vvqq_qq_m", "vvqq_qq_m", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="m_{qq} (GeV)", logY=True)
    makePlot("vvqq_qq_p", "vvqq_qq_p", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="p_{qq} (GeV)", logY=True)
    makePlot("vvqq_qq_costheta", "vvqq_qq_costheta", xMin=0, xMax=1, yMin=1e-1, yMax=-1, xLabel="vvqq_qq_costheta", logY=True)
    makePlot("vvqq_rp_no", "vvqq_rp_no", xMin=0, xMax=100, yMin=1e-1, yMax=-1, xLabel="vvqq_rp_no", logY=True)
    
    
    makePlot("vvqq_jet1_p", "vvqq_jet1_p", xMin=0, xMax=100, yMin=1e-1, yMax=-1, xLabel="vvqq_jet1_p", logY=True)
    makePlot("vvqq_jet2_p", "vvqq_jet2_p", xMin=0, xMax=100, yMin=1e-1, yMax=-1, xLabel="vvqq_jet2_p", logY=True)
    makePlot("vvqq_jet1_nc", "vvqq_jet1_nc", xMin=0, xMax=100, yMin=1e-1, yMax=-1, xLabel="vvqq_jet1_nc", logY=True)
    makePlot("vvqq_jet2_nc", "vvqq_jet2_nc", xMin=0, xMax=100, yMin=1e-1, yMax=-1, xLabel="vvqq_jet2_nc", logY=True)


    makePlot("vvqq_za_m", "vvqq_za_m", xMin=75, xMax=175, yMin=1e-1, yMax=-1, xLabel="vvqq_za_m", yLabel="Events", logY=True)
    makePlot("vvqq_za_p", "vvqq_za_p", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="vvqq_za_p", yLabel="Events", logY=True)
    makePlot("vvqq_recoil_za_m", "vvqq_recoil_za_m", xMin=75, xMax=175, yMin=1e-1, yMax=-1, xLabel="vv_qq_recoil_za_m", yLabel="Events", logY=True)
    makePlot("vvqq_recoil_z_m", "vvqq_recoil_z_m", xMin=75, xMax=175, yMin=1e-1, yMax=-1, xLabel="vv_qq_recoil_z_m", yLabel="Events", logY=True)

    makePlot("vvqq_vva_m", "vvqq_vva_m", xMin=75, xMax=175, yMin=1e-1, yMax=-1, xLabel="vvqq_vva_m", yLabel="Events", logY=True)
    makePlot("vvqq_vva_p", "vvqq_vva_p", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="vvqq_vva_p", yLabel="Events", logY=True)
    makePlot("vvqq_recoil_vva_m", "vvqq_recoil_vva_m", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="vvqq_recoil_vva_m", yLabel="Events", logY=True)
    makePlot("vvqq_recoil_vv_m", "vvqq_recoil_vv_m", xMin=75, xMax=175, yMin=1e-1, yMax=-1, xLabel="vvqq_recoil_vv_m", yLabel="Events", logY=True)

    

    makePlot("vvqq_za_p", "vvqq_za_p", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="vv_qq_za_p", yLabel="Events", logY=True)

    makePlot("vvqq_chi2", "vvqq_chi2", xMin=-100, xMax=100, yMin=1e-1, yMax=-1, xLabel="vv_qq_chi2", yLabel="Events", logY=True, rebin=10)
    makePlot("vvqq_mass_difference_qq", "vvqq_mass_difference_qq", xMin=0, xMax=50, yMin=1e-1, yMax=-1, xLabel="vvqq_mass_difference_qq", yLabel="Events", logY=True, rebin=2)
    makePlot("vvqq_mass_difference_vv", "vvqq_mass_difference_vv", xMin=0, xMax=50, yMin=1e-1, yMax=-1, xLabel="vvqq_mass_difference_vv", yLabel="Events", logY=True, rebin=2)



    outDir = "/home/submit/jaeyserm/public_html/fccee/h_za/plots/vvqq/qqa"
    makePlot("vvqqqqa_za_m", "vvqqqqa_za_m", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="vvqqqqa_za_m", logY=True)
    makePlot("vvqqqqa_recoil_za_m", "vvqqqqa_recoil_za_m", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="vvqqqqa_recoil_za_m", logY=True)
    makePlot("vvqqqqa_za_p", "vvqqqqa_za_p", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="vvqqqqa_za_p", logY=True)
    makePlot("vvqqqqa_final", "vvqqqqa_final", xMin=110, xMax=140, yMin=1e-1, yMax=-1, xLabel="vvqqqqa_final", logY=True)


    outDir = "/home/submit/jaeyserm/public_html/fccee/h_za/plots/vvqq/vva"
    makePlot("vvqqvva_z_p", "vvqqvva_z_p", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="vvqqvva_z_p", logY=True)
    makePlot("vvqqvva_recoil_z_m", "vvqqvva_recoil_z_m", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="vvqqvva_recoil_z_m", logY=True)
    makePlot("vvqqvva_final", "vvqqvva_final", xMin=110, xMax=140, yMin=1e-1, yMax=-1, xLabel="vvqqvva_final", logY=True)

    quit()

    #makePlot("tri_system_m_nOne", "tri_system_m_nOne", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="m_{qq#gamma} (GeV)", yLabel="Events", logY=True)
    #makePlot("tri_system_p_nOne", "tri_system_p_nOne", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="p_{qq#gamma} (GeV)", yLabel="Events", logY=True)
    
    
    #makePlot("recoil_qq_m", "recoil_qq_m", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="recoil qq", yLabel="Events", logY=True)
    #makePlot("recoil_qqa_m", "recoil_qqa_m", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="recoil qqa", yLabel="Events", logY=True)

    
    
    # vv_vv
    outDir = "/home/submit/jaeyserm/public_html/fccee/h_za/plots/vv_vv"
    makePlot("vv_vv_cosThetaMiss_nOne", "cosThetaMiss_nOne", xMin=0.9, xMax=1, yMin=1e-1, yMax=-1, xLabel="cos#theta_{miss}", yLabel="Events", logY=True, rebin=1)
    makePlot("vv_vv_missingEnergy_nOne", "missingEnergy_nOne", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="Missing energy (GeV)", yLabel="Events", logY=True, rebin=1)
    

    # vv_mumu
    outDir = "/home/submit/jaeyserm/public_html/fccee/h_za/plots/vv_mumu"
    makePlot("vv_mumu_mumu_m_nOne", "mumu_m_nOne", xMin=0, xMax=150, yMin=1e-1, yMax=-1, xLabel="m_{#mu#mu} (GeV)", yLabel="Events", logY=True)

    #makePlot("vv_qq_qq_m_nOne", "vv_qq_qq_m_nOne", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="m_{qq} (GeV)", yLabel="Events", logY=True)
    #makePlot("vv_qq_qq_p_nOne", "vv_qq_qq_p_nOne", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="p_{qq} (GeV)", yLabel="Events", logY=True)
    
    
    makePlot("vv_mumu_chi2", "vv_mumu_chi2", xMin=0, xMax=300, yMin=1e-1, yMax=-1, xLabel="vv_mumu_chi2", yLabel="Events", logY=True)
    makePlot("vv_mumu_za_m", "vv_mumu_za_m", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="vv_mumu_za_m", yLabel="Events", logY=True)
    makePlot("vv_mumu_recoil_za_m", "vv_mumu_recoil_za_m", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="vv_mumu_recoil_za_m", yLabel="Events", logY=True)
    makePlot("vv_mumu_recoil_z_m", "vv_mumu_recoil_z_m", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="vv_mumu_recoil_z_m", yLabel="Events", logY=True)





    # vv_ee
    outDir = "/home/submit/jaeyserm/public_html/fccee/h_za/plots/vv_ee"
    makePlot("vv_ee_ee_m_nOne", "ee_m_nOne", xMin=0, xMax=150, yMin=1e-1, yMax=-1, xLabel="m_{ee} (GeV)", yLabel="Events", logY=True)






    quit()
    makePlot("acoplanarity_nOne", "acoplanarity_nOne", xMin=0, xMax=2, yMin=1e-1, yMax=-1, xLabel="Acoplanarity", yLabel="Events", logY=True, rebin=1)
    makePlot("acolinearity_nOne", "acolinearity_nOne", xMin=0, xMax=2, yMin=1e-1, yMax=-1, xLabel="Acolinearity", yLabel="Events", logY=True, rebin=1)
    makePlot("acoplanarity", "acoplanarity", xMin=0, xMax=2, yMin=1e-1, yMax=-1, xLabel="Acoplanarity", yLabel="Events", logY=True, rebin=1)
    makePlot("acolinearity", "acolinearity", xMin=0, xMax=2, yMin=1e-1, yMax=-1, xLabel="Acolinearity", yLabel="Events", logY=True, rebin=1)
    makePlot("HZa_m", "HZa_m", xMin=110, xMax=130, yMin=0, yMax=-1, xLabel="m_{#gamma#gamma} (GeV)", yLabel="Events", logY=False)



    ### plots for each category
    makePlot("zqq_HZa_m", "zqq_HZa_m", xMin=110, xMax=130, yMin=0, yMax=-1, xLabel="m_{#gamma#gamma} (GeV)", yLabel="Events", logY=False)
    makePlot("zee_HZa_m", "zee_HZa_m", xMin=110, xMax=130, yMin=0, yMax=-1, xLabel="m_{#gamma#gamma} (GeV)", yLabel="Events", logY=False)
    makePlot("zmumu_HZa_m", "zmumu_HZa_m", xMin=110, xMax=130, yMin=0, yMax=-1, xLabel="m_{#gamma#gamma} (GeV)", yLabel="Events", logY=False)
    makePlot("znunu_HZa_m", "znunu_HZa_m", xMin=110, xMax=130, yMin=0, yMax=-1, xLabel="m_{#gamma#gamma} (GeV)", yLabel="Events", logY=False)

    quit()

    #
    makePlot("zqq_HZa_m", "zqq_HZa_m", xMin=110, xMax=130, yMin=0, yMax=1000, xLabel="m_{#gamma#gamma} (GeV)", yLabel="Events", logY=False, rebin=1)
    makePlot("zee_HZa_m", "zee_HZa_m", xMin=110, xMax=130, yMin=0, yMax=10, xLabel="m_{#gamma#gamma} (GeV)", yLabel="Events", logY=False, rebin=1)
    makePlot("zmumu_HZa_m", "zmumu_HZa_m", xMin=110, xMax=130, yMin=0, yMax=10, xLabel="m_{#gamma#gamma} (GeV)", yLabel="Events", logY=False, rebin=1)
    makePlot("znunu_HZa_m", "znunu_HZa_m", xMin=110, xMax=130, yMin=0, yMax=100, xLabel="m_{#gamma#gamma} (GeV)", yLabel="Events", logY=False, rebin=1)




    # N-1 plots
    #makePlot("zqq_HZa_m_nOne", "zqq_HZa_m_nOne", xMin=110, xMax=130, yMin=0, yMax=1000, xLabel="m_{#gamma#gamma} (GeV)", yLabel="Events", logY=False, rebin=1)




    

    
    
    makePlot("zqq_photons_p", "zqq_photons_p", xMin=0, xMax=150, yMin=0, yMax=2000, xLabel="m_{#gamma#gamma} (GeV)", yLabel="Events", logY=False, rebin=1)
    makePlot("zqq_photon_leading_p", "zqq_photon_leading_p", xMin=0, xMax=150, yMin=0, yMax=2000, xLabel="m_{#gamma#gamma} (GeV)", yLabel="Events", logY=False, rebin=1)
    makePlot("zqq_photons_p", "zqq_photons_p", xMin=0, xMax=150, yMin=0, yMax=2000, xLabel="m_{#gamma#gamma} (GeV)", yLabel="Events", logY=False, rebin=1)
    makePlot("zqq_photon_leading_costheta", "zqq_photon_leading_costheta", xMin=0, xMax=1, yMin=0, yMax=3000, xLabel="Acoplanarity", yLabel="Events", logY=False, rebin=50)
    makePlot("zqq_photons_costheta", "zqq_photons_costheta", xMin=0, xMax=1, yMin=0, yMax=3000, xLabel="Acoplanarity", yLabel="Events", logY=False, rebin=50)

    makePlot("zqq_acoplanarity", "zqq_acoplanarity", xMin=0, xMax=1, yMin=0, yMax=4000, xLabel="Acoplanarity", yLabel="Events", logY=False, rebin=10)
    makePlot("zqq_acolinearity", "zqq_acolinearity", xMin=0, xMax=1, yMin=0, yMax=4000, xLabel="Acolinearity", yLabel="Events", logY=False, rebin=10)
    
    makePlot("zqq_m_nOne", "zqq_m_nOne", xMin=50, xMax=120, yMin=0, yMax=5000, xLabel="m_{#gamma#gamma} (GeV)", yLabel="Events", logY=False, rebin=1)

    
    makePlot("photons_p", "photons_p", xMin=0, xMax=150, yMin=0, yMax=2000, xLabel="m_{#gamma#gamma} (GeV)", yLabel="Events", logY=False, rebin=1)
    
    # zqq plots\\
    makePlot("zqq_aa_p", "zqq_aa_p", xMin=0, xMax=100, yMin=0, yMax=2000, xLabel="m_{#gamma#gamma} (GeV)", yLabel="Events", logY=False, rebin=1)
    makePlot("zqq_qq_p", "zqq_qq_p", xMin=0, xMax=100, yMin=0, yMax=2000, xLabel="m_{#gamma#gamma} (GeV)", yLabel="Events", logY=False, rebin=1)
    
    makePlot("zqq_m", "zqq_m", xMin=0, xMax=120, yMin=0, yMax=2000, xLabel="m_{#gamma#gamma} (GeV)", yLabel="Events", logY=False, rebin=1)
    
    makePlot("zqq_acoplanarity", "zqq_acoplanarity", xMin=0, xMax=1, yMin=0, yMax=4000, xLabel="Acoplanarity", yLabel="Events", logY=False, rebin=10)
    makePlot("zqq_acolinearity", "zqq_acolinearity", xMin=0, xMax=1, yMin=0, yMax=4000, xLabel="Acolinearity", yLabel="Events", logY=False, rebin=10)
    
    
    
    quit()
    makePlot("mumu_recoil_m_nOne", "mumu_recoil_m_nOne", xMin=0, xMax=150, yMin=0, yMax=10000, xLabel="m_{rec} (GeV)", yLabel="Events", logY=False, rebin=1)
    makePlot("mumu_p_nOne", "mumu_p_nOne", xMin=0, xMax=100, yMin=0, yMax=10000, xLabel="p_{ll} (GeV)", yLabel="Events", logY=False, rebin=1)
    makePlot("cosThetaMiss_nOne", "cosThetaMiss_nOne", xMin=0, xMax=1, yMin=0, yMax=10000, xLabel="|cos(#theta_{miss})|", yLabel="Events", logY=False, rebin=100)
    makePlot("missingEnergy_nOne", "missingEnergy_nOne", xMin=0, xMax=150, yMin=0, yMax=10000, xLabel="Missing energy (GeV)", yLabel="Events", logY=False, rebin=1)