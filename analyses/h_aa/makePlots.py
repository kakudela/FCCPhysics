import sys
import os
import math
import copy
import array
import argparse

import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

sys.path.insert(0, f'{os.path.dirname(os.path.realpath(__file__))}/../../python')
import plotter


parser = argparse.ArgumentParser()
parser.add_argument("--run", help="Run combine", action='store_true')
parser.add_argument("--ecm", type=int, help="Center-of-mass energy", default=240)
args = parser.parse_args()



def getHist(hName, procs, rebin=1):
    hist = None
    for proc in procs:
        isVBF = False
        if 'nuenueVBF' in proc:
            proc = proc.replace('VBF', '')
            isVBF = True
        fIn = ROOT.TFile(f"{inputDir}/{proc}.root")
        h = fIn.Get(hName)
        h.SetDirectory(0)
        fIn.Close()
        if 'nuenue' in proc and isVBF:
            h_numunumu = getHist(hName, [proc.replace('nuenue', 'numunumu')], rebin=rebin)
            h.Add(h_numunumu, -1)
        if hist == None:
            hist = h
        else:
            hist.Add(h)
        
    hist.Rebin(rebin)
    #hist.Scale(3./10.8 * 1.05/1.39)
    return hist

def makeCutFlow(hName="cutFlow", cuts=[], labels=[], sig_scale=1.0):

    totEntries = 1 + len(bkgs)
    #leg = ROOT.TLegend(.5, 1.0-totEntries*0.06, .92, .90)
    leg = ROOT.TLegend(.45, 0.99-(len(bkgs)+2)*0.055, .95, .90)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.03)
    leg.SetMargin(0.2)

    hists_yields = []
    significances = []
    h_sig = getHist(hName, sigs)

    hists_yields.append(copy.deepcopy(h_sig))
    h_sig.Scale(sig_scale)
    h_sig.SetLineColor(sig_color)
    h_sig.SetLineWidth(4)
    h_sig.SetLineStyle(1)
    if sig_scale != 1:
        leg.AddEntry(h_sig, f"{sig_legend} (#times {int(sig_scale)})", "L")
    else:
        leg.AddEntry(h_sig, sig_legend, "L")

    # Get all bkg histograms
    st = ROOT.THStack()
    st.SetName("stack")
    h_bkg_tot = None
    for i,bkg in enumerate(bkgs):
        h_bkg = getHist(hName, bgks_cfg[bkg])

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

    for i,cut in enumerate(cuts):
        nsig = h_sig.GetBinContent(i+1)/sig_scale
        nbkg = 0
        for j,histProc in enumerate(hists_yields):
            nbkg = nbkg + histProc.GetBinContent(i+1)
            print(histProc.GetBinContent(i+1))
        print(nsig)
        s = nsig / (nsig + nbkg)**0.5
        print(i, cut, s)
        significances.append(s)

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

        'topRight'          : f"#sqrt{{s}} = {ecm} GeV, {lumi} ab^{{#minus1}}",
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
    canvas.SaveAs(f"{outDir}/{hName}.png")
    canvas.SaveAs(f"{outDir}/{hName}.pdf")

    out_orig = sys.stdout
    with open(f"{outDir}{hName}.txt", 'w') as f:
        sys.stdout = f

        formatted_row = '{:<10} {:<25} {:<25} {:<25} {:<25}' # adapt to #bkgs
        print(formatted_row.format(*(["Cut", "Significance", "Signal"]+bkgs)))
        print(formatted_row.format(*(["----------"]+["-----------------------"]*5)))
        for i,cut in enumerate(cuts):
            row = ["Cut %d"%i, "%.3f"%significances[i]]
            for j,histProc in enumerate(hists_yields):
                yield_, err = histProc.GetBinContent(i+1), histProc.GetBinError(i+1)
                row.append("%.2e +/- %.2e" % (yield_, err))

            print(formatted_row.format(*row))
    sys.stdout = out_orig

def makePlot(hName, outName, xMin=0, xMax=100, yMin=1, yMax=1e5, xLabel="xlabel", yLabel="Events", logX=False, logY=True, rebin=1, legPos=[0.3, 0.75, 0.9, 0.9], sig_scale=1):


    st = ROOT.THStack()
    st.SetName("stack")

    leg = ROOT.TLegend(legPos[0], legPos[1], legPos[2], legPos[3])
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.03)
    leg.SetMargin(0.2)

    h_sig = getHist(hName, sigs, rebin)
    if "TH2" in h_sig.ClassName(): h_sig = h_sig.ProjectionX("h_sig")
    h_sig.SetLineColor(sig_color)
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


def significance(hName, xMin=-10000, xMax=10000, reverse=False, sig_scale=1.):

    h_sig = getHist(hName, sigs)
    h_sig.Scale(sig_scale)
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
        if (sig+bkg) <= 0 or sig_tot <= 0:
            significance = 0
            sig_loss = 0
        else:
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

    if args.ecm == 240:
        ecm, lumi = 240, "10.8"
    if args.ecm == 365:
        ecm, lumi = 365, "3"

    inputDir = f"output/h_aa/histmaker/ecm{ecm}/"
    outDir = f"/home/submit/jaeyserm/public_html/fccee/h_aa/plots_ecm{ecm}/"


    bkgs = ["gaga", "Zgamma"]
    bkgs_legends = ["#gamma#gamma", "Z/#gamma^{*} #rightarrow f#bar{f}+#gamma(#gamma)"]
    # colors from https://github.com/mpetroff/accessible-color-cycles
    bkgs_colors = [ROOT.TColor.GetColor("#f89c20"), ROOT.TColor.GetColor("#5790fc")]
    sig_color = ROOT.TColor.GetColor("#e42536")

    if ecm == 240:
        sig_legend = "ZH(#gamma#gamma)"

        sigs = ['wzp6_ee_nunuH_Haa_ecm240', 'wzp6_ee_eeH_Haa_ecm240', 'wzp6_ee_tautauH_Haa_ecm240', 'wzp6_ee_ccH_Haa_ecm240', 'wzp6_ee_bbH_Haa_ecm240', 'wzp6_ee_qqH_Haa_ecm240', 'wzp6_ee_ssH_Haa_ecm240', 'wzp6_ee_mumuH_Haa_ecm240']
        #sigs = ['wzp6_ee_nunuH_Haa_ecm240']
        #sigs = ['wzp6_ee_ccH_Haa_ecm240', 'wzp6_ee_bbH_Haa_ecm240', 'wzp6_ee_qqH_Haa_ecm240', 'wzp6_ee_ssH_Haa_ecm240']
        #sigs = ['wzp6_ee_tautauH_Haa_ecm240']

        bgks_cfg = { 
            "gaga"      : ['wz3p6_ee_gammagamma_ecm240'],
            "Zgamma"    : ['wz3p6_ee_uu_ecm240', 'wz3p6_ee_dd_ecm240', 'wz3p6_ee_cc_ecm240', 'wz3p6_ee_ss_ecm240', 'wz3p6_ee_bb_ecm240', 'wz3p6_ee_tautau_ecm240', 'wz3p6_ee_mumu_ecm240', 'wz3p6_ee_nunu_ecm240', 'wz3p6_ee_ee_Mee_30_150_ecm240']
        }

        # cutflow
        cuts = ["cut0", "cut1", "cut2", "cut3", "cut4", "cut5", "cut6", "cut7", "cut8"]
        labels = ["All events", "#geq 2 #gamma", "#gamma ID", "80 < m_{rec} < 110", "20 < p_{#gamma#gamma} < 55", "0.2 < acol < 0.8", "acop > 0.05", "cos#theta_{miss} < 0.998", "120 < m_{#gamma#gamma} < 130"]
        makeCutFlow("cutFlow_zh", cuts, labels, 100.)


    if ecm == 365:
        sig_legend = "ZH(#gamma#gamma)"

        sigs = ['wzp6_ee_nunuH_Haa_ecm365', 'wzp6_ee_eeH_Haa_ecm365', 'wzp6_ee_tautauH_Haa_ecm365', 'wzp6_ee_ccH_Haa_ecm365', 'wzp6_ee_bbH_Haa_ecm365', 'wzp6_ee_qqH_Haa_ecm365', 'wzp6_ee_ssH_Haa_ecm365', 'wzp6_ee_mumuH_Haa_ecm365']
        #sigs = ['wzp6_ee_nunuH_Haa_ecm365']
        #sigs = ['wzp6_ee_ccH_Haa_ecm365', 'wzp6_ee_bbH_Haa_ecm365', 'wzp6_ee_qqH_Haa_ecm365', 'wzp6_ee_ssH_Haa_ecm365']
        #sigs = ['wzp6_ee_tautauH_Haa_ecm365']
        #sigs = ['wzp6_ee_nuenueVBFH_Haa_ecm365']
        #sigs = ['wzp6_ee_nuenueH_Haa_ecm365']
        #sigs = ['wzp6_ee_numunumuH_Haa_ecm365']
        
        bgks_cfg = { 
            "gaga"      : ['wz3p6_ee_gammagamma_ecm365'],
            "Zgamma"    : ['wz3p6_ee_uu_ecm365', 'wz3p6_ee_dd_ecm365', 'wz3p6_ee_cc_ecm365', 'wz3p6_ee_ss_ecm365', 'wz3p6_ee_bb_ecm365', 'wz3p6_ee_tautau_ecm365', 'wz3p6_ee_mumu_ecm365', 'wz3p6_ee_nunu_ecm365', 'wz3p6_ee_ee_Mee_30_150_ecm365']
        }

        # cutflow
        cuts = ["cut0", "cut1", "cut2", "cut3", "cut4", "cut5", "cut6", "cut7", "cut8"]
        labels = ["All events", "#geq 2 #gamma", "#gamma ID", "70 < m_{rec} < 230", "120 < p_{#gamma#gamma} < 150", "0.8 < acol < 1.75", "acop > 0.5", "cos#theta_{miss} < 0.999", "120 < m_{#gamma#gamma} < 130"]

        makeCutFlow("cutFlow_zh", cuts, labels, 100.)



    # significance
    if True and ecm == 240:
        significance("photon_leading_iso", 0, 0.5, reverse=True)
        significance("photon_leading_costheta", 0.5, 1, reverse=True)

        significance("photon_subleading_iso", 0, 2, reverse=True)
        significance("photon_subleading_costheta", 0.5, 1, reverse=True)
        
        significance("haa_recoil_m_nOne", 00, 100)
        significance("haa_recoil_m_nOne", 80, 300, reverse=True)
        
        significance("haa_p_nOne", 0, 100)
        significance("haa_p_nOne", 0, 100, reverse=True)
        
        significance("acolinearity_nOne", 0, 0.9)
        significance("acolinearity_nOne", 0, 0.9, reverse=True)
        
        significance("acoplanarity_nOne", 0, 1)
        significance("acoplanarity_nOne", 0, 1, reverse=True)
        
        significance("cosThetaMiss_nOne", 0.97, 1, reverse=True)
        ##significance("cosThetaPhotons", 0.5, 1, reverse=True)
        
        
        ## zqq optimization 
        significance("zqq_qq_p_nOne", 10, 70)
        significance("zqq_qq_p_nOne", 10, 70, reverse=True)

        significance("zqq_m_nOne", 0, 100)
        significance("zqq_m_nOne", 0, 120, reverse=True)

    # significance
    if True and ecm == 365:
        significance("photon_leading_iso", 0, 0.5, reverse=True)
        significance("photon_leading_costheta", 0.5, 1, reverse=True)

        significance("photon_subleading_iso", 0, 2, reverse=True)
        significance("photon_subleading_costheta", 0.5, 1, reverse=True)
        
        significance("haa_recoil_m_nOne", 00, 100)
        significance("haa_recoil_m_nOne", 80, 300, reverse=True)
        
        significance("haa_p_nOne", 0, 150)
        significance("haa_p_nOne", 50, 170, reverse=True)
        
        significance("acolinearity_nOne", 0, 2)
        significance("acolinearity_nOne", 0, 2, reverse=True)
        
        significance("acoplanarity_nOne", 0, 3)
        significance("acoplanarity_nOne", 0, 3, reverse=True)
        
        significance("cosThetaMiss_nOne", 0.97, 1, reverse=True)
        ##significance("cosThetaPhotons", 0.5, 1, reverse=True)
        
        
        ## zqq optimization 
        significance("zqq_qq_p_nOne", 20, 160)
        significance("zqq_qq_p_nOne", 20, 120, reverse=True)

        significance("zqq_m_nOne", 0, 100)
        significance("zqq_m_nOne", 80, 110, reverse=True)
        
        



    makePlot("photons_all_p", "photons_all_p", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="Photons momentum (GeV)", yLabel="Events", logY=True)
    makePlot("photons_p", "photons_p", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="Photons momentum (GeV)", yLabel="Events", logY=True)
    makePlot("photons_n", "photons_n", xMin=0, xMax=15, yMin=1e-1, yMax=-1, xLabel="Number of photons", yLabel="Events", logY=True)
    makePlot("photon_leading_p", "photon_leading_p", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="Leading photon momentum (GeV)", yLabel="Events", logY=True)
    makePlot("photon_subleading_p", "photon_subleading_p", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="Subleading photon momentum (GeV)", yLabel="Events", logY=True)

    makePlot("photon_leading_iso", "photon_leading_iso", xMin=0, xMax=2, yMin=1e-1, yMax=-1, xLabel="Leading photon isolation", yLabel="Events", logY=True)
    makePlot("photon_subleading_iso", "photon_subleading_iso", xMin=0, xMax=2, yMin=1e-1, yMax=-1, xLabel="Subleading photon isolation", yLabel="Events", logY=True)


    makePlot("photon_leading_costheta", "photon_leading_costheta", xMin=0, xMax=1, yMin=1e-1, yMax=-1, xLabel="cos#theta leading photon", yLabel="Events", logY=True)
    makePlot("photon_subleading_costheta", "photon_subleading_costheta", xMin=0, xMax=1, yMin=1e-1, yMax=-1, xLabel="cos#theta subleading photon", yLabel="Events", logY=True)


    makePlot("photons_n", "photons_n", xMin=0, xMax=10, yMin=1e-1, yMax=-1, xLabel="Photon multiplicity", yLabel="Events", logY=True)

    makePlot("haa_recoil_m_nOne", "haa_recoil_m_nOne", xMin=0, xMax=400, yMin=1e-1, yMax=-1, xLabel="Recoil", yLabel="Events", logY=True, rebin=1)
    makePlot("haa_p_nOne", "haa_p_nOne", xMin=0, xMax=400, yMin=1e-1, yMax=-1, xLabel="p_{#gamma#gamma} (GeV)", yLabel="Events", logY=True, rebin=1)
    makePlot("acoplanarity_nOne", "acoplanarity_nOne", xMin=0, xMax=5, yMin=1e-1, yMax=-1, xLabel="Acoplanarity", yLabel="Events", logY=True, rebin=1)
    makePlot("acolinearity_nOne", "acolinearity_nOne", xMin=0, xMax=2, yMin=1e-1, yMax=-1, xLabel="Acolinearity", yLabel="Events", logY=True, rebin=1)
    makePlot("cosThetaMiss_nOne", "cosThetaMiss_nOne", xMin=0.9, xMax=1, yMin=1e-1, yMax=-1, xLabel="cos#theta_{miss}", yLabel="Events", logY=True, rebin=1)
    makePlot("missingEnergy", "missingEnergy", xMin=0, xMax=150, yMin=1e-1, yMax=-1, xLabel="Missing energy (GeV)", yLabel="Events", logY=True, rebin=1)
    makePlot("missingMass", "missingMass", xMin=0, xMax=300, yMin=1e-1, yMax=-1, xLabel="Missing mass (GeV)", yLabel="Events", logY=True, rebin=1)
    makePlot("haa_m_nOne", "haa_m_nOne", xMin=120, xMax=130, yMin=1e-1, yMax=-1, xLabel="m_{#gamma#gamma} (GeV)", yLabel="Events", logY=True, rebin=1)
    makePlot("haa_m_nOne", "haa_m_nOne_nolog", xMin=120, xMax=130, yMin=0, yMax=-1, xLabel="m_{#gamma#gamma} (GeV)", yLabel="Events", logY=False, rebin=1)


    makePlot("photon_momentum_ratio", "photon_momentum_ratio", xMin=0, xMax=5, yMin=1e-1, yMax=-1, xLabel="Energy ratio", yLabel="Events", logY=True, rebin=1)

    makePlot("photon_dphi", "photon_dphi", xMin=-4, xMax=4, yMin=1e-1, yMax=-1, xLabel="photon_dphi", yLabel="Events", logY=True, rebin=1)
    makePlot("cosThetaPhotons", "cosThetaPhotons", xMin=0, xMax=1, yMin=1e-1, yMax=-1, xLabel="cosThetaPhotons_nOne", yLabel="Events", logY=True, rebin=1)

    makePlot("haa_costheta", "haa_costheta", xMin=0, xMax=1, yMin=1e-1, yMax=-1, xLabel="haa_costheta", yLabel="Events", logY=True, rebin=1)





    # mumu
    outDir_base = outDir
    outDir = f"{outDir_base}/mumu/"
    makePlot("muons_all_p", "muons_all_p", xMin=0, xMax=100, yMin=1e-1, yMax=-1, xLabel="Muon momentum (GeV)", yLabel="Events", logY=True, rebin=1)
    makePlot("zmumu_m_nOne", "zmumu_m_nOne", xMin=50, xMax=120, yMin=1e-1, yMax=-1, xLabel="Muon mass (GeV)", yLabel="Events", logY=True, rebin=1)
    makePlot("zmumu_p_nOne", "zmumu_p_nOne", xMin=0, xMax=120, yMin=1e-1, yMax=-1, xLabel="Muon p (GeV)", yLabel="Events", logY=True, rebin=1)
    makePlot("zmumu_haa_m", "zmumu_haa_m", xMin=120, xMax=130, yMin=1e-1, yMax=-1, xLabel="m_{#gamma#gamma} (GeV)", yLabel="Events", logY=True, rebin=1)
    makePlot("zmumu_haa_m", "zmumu_haa_m_nolog", xMin=120, xMax=130, yMin=0, yMax=-1, xLabel="m_{#gamma#gamma} (GeV)", yLabel="Events", logY=False, rebin=1)
    makePlot("zmumu_dtheta", "zmumu_dtheta", xMin=-4, xMax=4, yMin=1e-1, yMax=-1, xLabel="dtheta", yLabel="Events", logY=True, rebin=1)
    makePlot("zmumu_dphi", "zmumu_dphi", xMin=-4, xMax=4, yMin=1e-1, yMax=-1, xLabel="dphi", yLabel="Events", logY=True, rebin=1)

    makePlot("zmumu_costheta1", "zmumu_costheta1", xMin=-1, xMax=1, yMin=1e-1, yMax=-1, xLabel="zmumu_costheta1", yLabel="Events", logY=True, rebin=1)
    makePlot("zmumu_costheta2", "zmumu_costheta2", xMin=-1, xMax=1, yMin=1e-1, yMax=-1, xLabel="zmumu_costheta2", yLabel="Events", logY=True, rebin=1)
    makePlot("zmumu_photon_mumu_costheta", "zmumu_photon_mumu_costheta", xMin=-1, xMax=1, yMin=1e-1, yMax=-1, xLabel="zmumu_photon_mumu_costheta", yLabel="Events", logY=True, rebin=1)

    # ee
    outDir = f"{outDir_base}/ee/"
    makePlot("electrons_all_p", "electrons_all_p", xMin=0, xMax=100, yMin=1e-1, yMax=-1, xLabel="Electron momentum (GeV)", yLabel="Events", logY=True, rebin=1)
    makePlot("zee_m_nOne", "zee_m_nOne", xMin=50, xMax=120, yMin=1e-1, yMax=-1, xLabel="Electron mass (GeV)", yLabel="Events", logY=True, rebin=1)
    #makePlot("zee_p_nOne", "zee_p_nOne", xMin=0, xMax=120, yMin=1e-1, yMax=-1, xLabel="Electron p (GeV)", yLabel="Events", logY=True, rebin=1)
    makePlot("zee_haa_m", "zee_haa_m", xMin=120, xMax=130, yMin=1e-1, yMax=-1, xLabel="m_{#gamma#gamma} (GeV)", yLabel="Events", logY=True, rebin=1)
    makePlot("zee_haa_m", "zee_haa_m_nolog", xMin=120, xMax=130, yMin=0, yMax=-1, xLabel="m_{#gamma#gamma} (GeV)", yLabel="Events", logY=False, rebin=1)

    # qq
    outDir = f"{outDir_base}/qq/" 

    makePlot("zqq_jet1_p_nOne", "zqq_jet1_p_nOne", xMin=0, xMax=200, yMin=1e-5, yMax=-1, xLabel="zqq_jet1_p_nOne (GeV)", yLabel="Events", logY=True, rebin=1)
    makePlot("zqq_jet2_p_nOne", "zqq_jet2_p_nOne", xMin=0, xMax=200, yMin=1e-5, yMax=-1, xLabel="zqq_jet2_p_nOne (GeV)", yLabel="Events", logY=True, rebin=1)
    makePlot("zqq_qq_p_nOne", "zqq_qq_p_nOne", xMin=0, xMax=200, yMin=1e-5, yMax=-1, xLabel="p_qq (GeV)", yLabel="Events", logY=True, rebin=1)
    makePlot("zqq_m_nOne", "zqq_m_nOne", xMin=50, xMax=150, yMin=1e-1, yMax=-1, xLabel="z_qq (GeV)", yLabel="Events", logY=True, rebin=1)
    makePlot("zqq_haa_m", "zqq_haa_m", xMin=120, xMax=130, yMin=1e-1, yMax=-1, xLabel="m_{#gamma#gamma} (GeV)", yLabel="Events", logY=True, rebin=1)
    makePlot("zqq_haa_m", "zqq_haa_m_nolog", xMin=120, xMax=130, yMin=0, yMax=-1, xLabel="m_{#gamma#gamma} (GeV)", yLabel="Events", logY=False, rebin=1)
    
    makePlot("zqq_costheta1", "zqq_costheta1", xMin=0, xMax=1, yMin=1e-1, yMax=-1, xLabel="zqq_costheta1", yLabel="Events", logY=True, rebin=1)
    makePlot("zqq_costheta2_1", "zqq_costheta2_1", xMin=0, xMax=1, yMin=1e-1, yMax=-1, xLabel="zqq_costheta2_1", yLabel="Events", logY=True, rebin=1)
    makePlot("zqq_costheta2_2", "zqq_costheta2_2", xMin=0, xMax=1, yMin=1e-1, yMax=-1, xLabel="zqq_costheta2_2", yLabel="Events", logY=True, rebin=1)
    makePlot("zqq_photon_qq_costheta", "zqq_photon_qq_costheta", xMin=0, xMax=1, yMin=1e-1, yMax=-1, xLabel="zqq_photon_qq_costheta", yLabel="Events", logY=True, rebin=1)



    # nunu
    outDir = f"{outDir_base}/nunu/"
    makePlot("znunu_haa_m", "znunu_haa_m", xMin=120, xMax=130, yMin=1e-1, yMax=-1, xLabel="m_{#gamma#gamma} (GeV)", yLabel="Events", logY=True, rebin=1)
    makePlot("znunu_haa_m", "znunu_haa_m_nolog", xMin=120, xMax=130, yMin=0, yMax=-1, xLabel="m_{#gamma#gamma} (GeV)", yLabel="Events", logY=False, rebin=1)
    makePlot("znunu_cosThetaMiss", "znunu_cosThetaMiss_nOne", xMin=0, xMax=1, yMin=1e-1, yMax=-1, xLabel="cos#theta_{miss}", yLabel="Events", logY=True, rebin=100)
    makePlot("znunu_missingEnergy", "znunu_missingEnergy_nOne", xMin=80, xMax=120, yMin=1e-1, yMax=-1, xLabel="Missing energy (GeV)", yLabel="Events", logY=True, rebin=1)
    makePlot("znunu_cosThetaPhotons", "znunu_cosThetaPhotons", xMin=0, xMax=1, yMin=1e-1, yMax=-1, xLabel="znunu_cosThetaPhotons", yLabel="Events", logY=True, rebin=1)


    # vbf
    if ecm == 365:

        outDir = f"{outDir_base}/vbf/"

        sig_legend = "#nu#bar{#nu}H(#gamma#gamma)"
        sigs = ['wzp6_ee_nuenueVBFH_Haa_ecm365']

        cuts_vbf = ["cut0", "cut1", "cut2", "cut3", "cut4", "cut5", "cut6", "cut7"]
        labels_vbf = ["All events", "#geq 2 #gamma", "#gamma ID", "m_{rec} < 240", "30 < p_{#gamma#gamma} > 130", "0.35 < acol > 1.6", "acop > 0.2", "Missing mass > 170", "120 < m_{#gamma#gamma} < 130"]
        makeCutFlow("cutFlow_vbf", cuts_vbf, labels_vbf, 100.)


        ## vbf 
        significance("vbf_haa_recoil_m_nOne", 0, 250)
        significance("vbf_haa_recoil_m_nOne", 0, 250, reverse=True)
        
        significance("vbf_haa_p_nOne", 0, 200, sig_scale=3.)
        significance("vbf_haa_p_nOne", 0, 200, reverse=True)
        
        significance("vbf_acolinearity_nOne", 0, 4, sig_scale=3.)
        significance("vbf_acolinearity_nOne", 0, 4, reverse=True)
        
        significance("vbf_acoplanarity_nOne", 0, 4)
        significance("vbf_acoplanarity_nOne", 0, 4, reverse=True)
        
        significance("vbf_haa_m_nOne", 100, 150)
        significance("vbf_haa_m_nOne", 100, 150, reverse=True)

        significance("vbf_cosThetaMiss_nOne", 0.95, 1, reverse=True)
        significance("vbf_missingMass_nOne", 100, 250)
        
        makePlot("vbf_haa_recoil_m_nOne", "vbf_haa_recoil_m_nOne", xMin=0, xMax=400, yMin=1e-1, yMax=-1, xLabel="Recoil", yLabel="Events", logY=True, rebin=1)
        makePlot("vbf_haa_p_nOne", "vbf_haa_p_nOne", xMin=0, xMax=400, yMin=1e-1, yMax=-1, xLabel="p_{#gamma#gamma} (GeV)", yLabel="Events", logY=True, rebin=1)
        makePlot("vbf_acoplanarity_nOne", "vbf_acoplanarity_nOne", xMin=0, xMax=5, yMin=1e-1, yMax=-1, xLabel="Acoplanarity", yLabel="Events", logY=True, rebin=1)
        makePlot("vbf_acolinearity_nOne", "vbf_acolinearity_nOne", xMin=0, xMax=2, yMin=1e-1, yMax=-1, xLabel="Acolinearity", yLabel="Events", logY=True, rebin=1)
        makePlot("vbf_cosThetaMiss_nOne", "vbf_cosThetaMiss_nOne", xMin=0, xMax=1, yMin=1e-1, yMax=-1, xLabel="cos#theta_{miss}", yLabel="Events", logY=True, rebin=100)
        makePlot("vbf_haa_m_nOne", "vbf_haa_m_nOne", xMin=120, xMax=130, yMin=1e-1, yMax=-1, xLabel="m_{#gamma#gamma} (GeV)", yLabel="Events", logY=True, rebin=1)
        makePlot("vbf_haa_m_nOne", "vbf_haa_m_nOne_nolog", xMin=120, xMax=130, yMin=0, yMax=-1, xLabel="m_{#gamma#gamma} (GeV)", yLabel="Events", logY=False, rebin=5)
        makePlot("vbf_missingMass_nOne", "vbf_missingMass_nOne", xMin=0, xMax=400, yMin=1e-1, yMax=-1, xLabel="p_{#gamma#gamma} (GeV)", yLabel="Events", logY=True, rebin=1)

        makePlot("vbf_haa_costheta_nOne", "vbf_haa_costheta_nOne", xMin=0, xMax=1, yMin=1e-1, yMax=-1, xLabel="vbf_haa_costheta_nOne", yLabel="Events", logY=True, rebin=1)

