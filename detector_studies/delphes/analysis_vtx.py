

import sys, os, glob, math
import ROOT
import logging

logging.basicConfig(format='%(levelname)s: %(message)s')
logger = logging.getLogger("fcclogger")
logger.setLevel(logging.INFO)


ROOT.gROOT.SetBatch()
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

ROOT.EnableImplicitMT() # use all cores
#ROOT.DisableImplicitMT() # single core

# load libraries
ROOT.gSystem.Load("libFCCAnalyses")
fcc_loaded = ROOT.dummyLoader()
ROOT.gInterpreter.Declare("using namespace FCCAnalyses;")
ROOT.gInterpreter.Declare('#include "functions.h"')



bins_p = (250, 0, 250)
bins_res = (10000, -0.1, 0.1)
bins_counts = (100, 0, 100)
bins_distance_abs = (10000, 0, 10)
bins_distance = (100000, -50, 50)

## functions defined here: https://github.com/HEP-FCC/FCCAnalyses/blob/master/analyzers/dataframe/src/myUtils.cc
## vertex fitter: https://indico.cern.ch/event/1003610/contributions/4214579/attachments/2187815/3696958/Bedeschi_Vertexing_Feb2021.pdf
## perf. plots: https://indico.cern.ch/event/965346/contributions/4062989/attachments/2125687/3578824/vertexing.pdf

def analysis(input_files, output_file):

    df = ROOT.RDataFrame("events", input_files)

    df = df.Alias("MCRecoAssociations0", "_MCRecoAssociations_rec.index")
    df = df.Alias("MCRecoAssociations1", "_MCRecoAssociations_sim.index")
    df = df.Alias("Particle0", "_Particle_parents.index")

    df = df.Define("MC_p", "FCCAnalyses::MCParticle::get_p(Particle)")
    df = df.Define("MC_theta", "FCCAnalyses::MCParticle::get_theta(Particle)")
    df = df.Define("MC_charge", "FCCAnalyses::MCParticle::get_charge(Particle)")
    df = df.Define("MC_phi", "FCCAnalyses::MCParticle::get_phi(Particle)")

    df = df.Define("RP_p", "ReconstructedParticle::get_p(ReconstructedParticles)")
    df = df.Define("RP_theta", "ReconstructedParticle::get_theta(ReconstructedParticles)")
    df = df.Define("RP_charge", "ReconstructedParticle::get_charge(ReconstructedParticles)")
    df = df.Define("RP_phi", "ReconstructedParticle::get_phi(ReconstructedParticles)")


    # MC Vertex
    df = df.Define("MCVertexObject", "myUtils::get_MCVertexObject(Particle, Particle0)")
    df = df.Define("MC_Vertex_x", "myUtils::get_MCVertex_x(MCVertexObject)")
    df = dfdf = df.Define("MC_Vertex_y", "myUtils::get_MCVertex_y(MCVertexObject)")
    df = df.Define("MC_Vertex_z", "myUtils::get_MCVertex_z(MCVertexObject)")
    df = df.Define("MC_Vertex_ind", "myUtils::get_MCindMCVertex(MCVertexObject)")
    df = df.Define("MC_Vertex_ntrk", "myUtils::get_NTracksMCVertex(MCVertexObject)")
    df = df.Define("MC_Vertex_n", "int(MC_Vertex_x.size())")




    # RECO Vertex
    df = df.Define("VertexObject", "myUtils::get_VertexObject(MCVertexObject, ReconstructedParticles, _EFlowTrack_trackStates, MCRecoAssociations0, MCRecoAssociations1)")
    df = df.Define("Vertex_x", "myUtils::get_Vertex_x(VertexObject)")
    df = df.Define("Vertex_y", "myUtils::get_Vertex_y(VertexObject)")
    df = df.Define("Vertex_z", "myUtils::get_Vertex_z(VertexObject)")
    df = df.Define("Vertex_xErr", "myUtils::get_Vertex_xErr(VertexObject)")
    df = df.Define("Vertex_yErr", "myUtils::get_Vertex_yErr(VertexObject)")
    df = df.Define("Vertex_zErr", "myUtils::get_Vertex_zErr(VertexObject)")

    df = df.Define("Vertex_chi2", "myUtils::get_Vertex_chi2(VertexObject)")
    #df = df.Define("Vertex_mcind", "myUtils::get_Vertex_indMC(VertexObject, MCVertexObject)")
    df = df.Define("Vertex_mcind", "myUtils::get_Vertex_indMC(VertexObject)")
    df = df.Define("Vertex_ind", "myUtils::get_Vertex_ind(VertexObject)")
    df = df.Define("Vertex_isPV", "myUtils::get_Vertex_isPV(VertexObject)") # is primary vertex PV?
    df = df.Define("Vertex_ntrk", "myUtils::get_Vertex_ntracks(VertexObject)") # number of tracks per vertex
    df = df.Define("Vertex_n", "int(Vertex_x.size())") # number of RECO vertices

    # distance of any vertex to the primary vertex (PV)
    df = df.Define("Vertex_d2PV",  "myUtils::get_Vertex_d2PV(VertexObject, -1)")
    df = df.Define("Vertex_d2PVx", "myUtils::get_Vertex_d2PV(VertexObject, 0)")
    df = df.Define("Vertex_d2PVy", "myUtils::get_Vertex_d2PV(VertexObject, 1)")
    df = df.Define("Vertex_d2PVz", "myUtils::get_Vertex_d2PV(VertexObject, 2)")

    df = df.Define("Vertex_d2PVErr", "myUtils::get_Vertex_d2PVError(VertexObject, -1)")
    df = df.Define("Vertex_d2PVxErr","myUtils::get_Vertex_d2PVError(VertexObject, 0)")
    df = df.Define("Vertex_d2PVyErr","myUtils::get_Vertex_d2PVError(VertexObject, 1)")
    df = df.Define("Vertex_d2PVzErr","myUtils::get_Vertex_d2PVError(VertexObject, 2)")

    df = df.Define("Vertex_d2PVSig",  "Vertex_d2PV/Vertex_d2PVErr")
    df = df.Define("Vertex_d2PVxSig", "Vertex_d2PVx/Vertex_d2PVxErr")
    df = df.Define("Vertex_d2PVySig", "Vertex_d2PVy/Vertex_d2PVyErr")
    df = df.Define("Vertex_d2PVzSig", "Vertex_d2PVz/Vertex_d2PVzErr")

    # distance of MC vs reco vertex (in um)
    df = df.Define("Vertex_d2MC", "myUtils::get_Vertex_d2MC(VertexObject, MCVertexObject, Vertex_mcind, -1)*1000")
    df = df.Define("Vertex_d2MCx", "myUtils::get_Vertex_d2MC(VertexObject, MCVertexObject, Vertex_mcind, 0)*1000")
    df = df.Define("Vertex_d2MCy", "myUtils::get_Vertex_d2MC(VertexObject, MCVertexObject, Vertex_mcind, 1)*1000")
    df = df.Define("Vertex_d2MCz", "myUtils::get_Vertex_d2MC(VertexObject, MCVertexObject, Vertex_mcind, 2)*1000")

    df = df.Define("VertexPV_d2MC", "Vertex_d2MC[Vertex_isPV]")
    df = df.Define("VertexPV_d2MCx", "Vertex_d2MCx[Vertex_isPV]")
    df = df.Define("VertexPV_d2MCy", "Vertex_d2MCy[Vertex_isPV]")
    df = df.Define("VertexPV_d2MCz", "Vertex_d2MCz[Vertex_isPV]")

    h_MC_Vertex_ntrk = df.Histo1D(("MC_Vertex_ntrk", "", *bins_counts), "MC_Vertex_ntrk")
    h_MC_Vertex_n = df.Histo1D(("MC_Vertex_n", "", *bins_counts), "MC_Vertex_n")

    h_Vertex_ntrk = df.Histo1D(("Vertex_ntrk", "", *bins_counts), "Vertex_ntrk")
    h_Vertex_n = df.Histo1D(("Vertex_n", "", *bins_counts), "Vertex_n")


    h_Vertex_d2MC = df.Histo1D(("Vertex_d2MC", "", *bins_distance_abs), "Vertex_d2MC")
    h_VertexPV_d2MC = df.Histo1D(("VertexPV_d2MC", "", *bins_distance_abs), "VertexPV_d2MC")
    h_VertexPV_d2MCx = df.Histo1D(("VertexPV_d2MCx", "", *bins_distance), "VertexPV_d2MCx")
    h_VertexPV_d2MCy = df.Histo1D(("VertexPV_d2MCy", "", *bins_distance), "VertexPV_d2MCy")
    h_VertexPV_d2MCz = df.Histo1D(("VertexPV_d2MCz", "", *bins_distance), "VertexPV_d2MCz")
    h_Vertex_d2PV = df.Histo1D(("Vertex_d2PV", "", *bins_distance_abs), "Vertex_d2PV")


    h_Vertex_isPV = df.Histo1D(("Vertex_isPV", "", *bins_counts), "Vertex_isPV")


    # write output
    fout = ROOT.TFile(output_file, "RECREATE")
    h_MC_Vertex_ntrk.Write()
    h_MC_Vertex_n.Write()
    h_Vertex_ntrk.Write()
    h_Vertex_n.Write()
    h_Vertex_d2MC.Write()
    h_VertexPV_d2MC.Write()
    h_VertexPV_d2MCx.Write()
    h_VertexPV_d2MCy.Write()
    h_VertexPV_d2MCz.Write()
    h_Vertex_d2PV.Write()
    h_Vertex_isPV.Write()
    fout.Close()



if __name__ == "__main__":

    input_files, output_file = ["samples/CLD_2T_Zbb_ecm91p2.root"], "output/vtx_CLD_2T_Zbb_ecm91p2.root"

    logger.info(f"Start analysis")
    analysis(input_files, output_file)
    logger.info(f"Done! Output saved to {output_file}")

