
import ROOT
import array
import numpy as np
import socket

ROOT.TH1.SetDefaultSumw2(ROOT.kTRUE)
from addons.TMVAHelper.TMVAHelper import TMVAHelperXGB


# this number is used to run only on a fraction of the Monte Carlo samples, which can speed up the analysis
fraction = 0.02

processList = {
    'p8_ee_WW_mumu_ecm240':             {'fraction':fraction},
    'p8_ee_ZZ_ecm240':                  {'fraction':fraction},
    'wzp6_ee_tautau_ecm240':            {'fraction':fraction},
    'wzp6_ee_mumu_ecm240':              {'fraction':fraction},

    'wzp6_ee_mumuH_Hbb_ecm240':         {'fraction':fraction},
    'wzp6_ee_mumuH_Hcc_ecm240':         {'fraction':fraction},
    'wzp6_ee_mumuH_Hss_ecm240':         {'fraction':fraction},
    'wzp6_ee_mumuH_Hgg_ecm240':         {'fraction':fraction},
    'wzp6_ee_mumuH_Haa_ecm240':         {'fraction':fraction},
    'wzp6_ee_mumuH_HZa_ecm240':         {'fraction':fraction},
    'wzp6_ee_mumuH_HWW_ecm240':         {'fraction':fraction},
    'wzp6_ee_mumuH_HZZ_ecm240':         {'fraction':fraction},
    'wzp6_ee_mumuH_Hmumu_ecm240':       {'fraction':fraction},
    'wzp6_ee_mumuH_Htautau_ecm240':     {'fraction':fraction},
}




if "mit.edu" in socket.gethostname(): # configuration for MIT
    inputDir = "/ceph/submit/data/group/fcc/ee/generation/DelphesEvents/winter2023/IDEA/"
    procDict = "/ceph/submit/data/group/fcc/ee/generation/DelphesEvents/winter2023/IDEA/samplesDict.json"
else: # default configuration for CERN
    prodTag     = "FCCee/winter2023/IDEA/"
    procDict = "FCCee_procDict_winter2023_IDEA.json"

# additional/custom C++ functions
includePaths = ["utils.h"]


# output directory
outputDir   = "output/ZmumuHbb/histmaker/"

# optional: ncpus, default is 4, -1 uses all cores available
nCPUS       = 8


# scale the histograms with the cross-section and integrated luminosity
doScale = True
intLumi = 10.8e6 # integrated luminosity at 240 GeV


# define histograms
bins_m = (250, 0, 250) # 100 MeV bins
bins_cosThetaMiss = (10000, 0, 1)

bins_theta = (400, 0, 4)
bins_phi = (400, -4, 4)

bins_count = (50, 0, 50)
bins_charge = (10, -5, 5)


################################################################################
## load modules and files for jet clustering and flavor tagging
################################################################################

if "mit.edu" in socket.gethostname(): # configuration for MIT
    weaver_preproc = "/ceph/submit/data/group/fcc/ee/aux/tagger/fccee_flavtagging_edm4hep_wc.json"
    weaver_model = "/ceph/submit/data/group/fcc/ee/aux/tagger/fccee_flavtagging_edm4hep_wc.onnx"
else: # default configuration for CERN
    weaver_preproc = "/eos/experiment/fcc/ee/jet_flavour_tagging/winter2023/wc_pt_7classes_12_04_2023/fccee_flavtagging_edm4hep_wc.json"
    weaver_model = "/eos/experiment/fcc/ee/jet_flavour_tagging/winter2023/wc_pt_7classes_12_04_2023/fccee_flavtagging_edm4hep_wc.onnx"

from addons.ONNXRuntime.jetFlavourHelper import JetFlavourHelper
from addons.FastJet.jetClusteringHelper import ExclusiveJetClusteringHelper
from examples.FCCee.weaver.config import collections, njets

################################################################################


def build_graph(df, dataset):

    hists = []

    df = df.Define("weight", "1.0")
    weightsum = df.Sum("weight")

    df = df.Define("cut0", "0")
    df = df.Define("cut1", "1")
    df = df.Define("cut2", "2")
    df = df.Define("cut3", "3")
    df = df.Define("cut4", "4")
    df = df.Define("cut5", "5")
    df = df.Define("cut6", "6")
    df = df.Define("cut7", "7")
    df = df.Define("cut8", "8")
    df = df.Define("cut9", "9")
    df = df.Define("cut10", "10")
    df = df.Define("cut11", "11")
    df = df.Define("cut12", "12")
    df = df.Define("cut13", "13")
    df = df.Define("cut14", "14")

    df = df.Alias("Particle0", "Particle#0.index")
    df = df.Alias("Particle1", "Particle#1.index")
    df = df.Alias("MCRecoAssociations0", "MCRecoAssociations#0.index")
    df = df.Alias("MCRecoAssociations1", "MCRecoAssociations#1.index")
    df = df.Alias("Muon0", "Muon#0.index")

    # Get all the reconstructed muons
    df = df.Define("muons_all", "FCCAnalyses::ReconstructedParticle::get(Muon0, ReconstructedParticles)")
    df = df.Define("muons_all_p", "FCCAnalyses::ReconstructedParticle::get_p(muons_all)")

    # Select muons that have at least 20 GeV in momentum
    df = df.Define("muons", "FCCAnalyses::ReconstructedParticle::sel_p(20)(muons_all)")
    df = df.Define("muons_p", "FCCAnalyses::ReconstructedParticle::get_p(muons)")
    df = df.Define("muons_theta", "FCCAnalyses::ReconstructedParticle::get_theta(muons)")
    df = df.Define("muons_phi", "FCCAnalyses::ReconstructedParticle::get_phi(muons)")
    df = df.Define("muons_q", "FCCAnalyses::ReconstructedParticle::get_charge(muons)")
    df = df.Define("muons_no", "FCCAnalyses::ReconstructedParticle::get_n(muons)")


    # baseline selections and histograms
    hists.append(df.Histo1D(("muons_all_p", "", *bins_m), "muons_all_p"))

    hists.append(df.Histo1D(("muons_p", "", *bins_m), "muons_p"))
    hists.append(df.Histo1D(("muons_theta", "", *bins_theta), "muons_theta"))
    hists.append(df.Histo1D(("muons_phi", "", *bins_phi), "muons_phi"))
    hists.append(df.Histo1D(("muons_q", "", *bins_charge), "muons_q"))
    hists.append(df.Histo1D(("muons_no", "", *bins_count), "muons_no"))


    #########
    ### CUT 0 all events
    #########
    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut0"))


    #########
    ### CUT 1: 2 opposite-sign muons
    #########
    df = df.Filter("muons_no == 2 && Sum(muons_q) == 0")
    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut1"))

    # Construct the Lorentz vectors of the muons
    df = df.Define("muons_tlv", "FCCAnalyses::ReconstructedParticle::get_tlv(muons)")
    df = df.Define("muon1_tlv", "muons_tlv[0]")
    df = df.Define("muon2_tlv", "muons_tlv[1]")
    df = df.Define("zmumu", "muon1_tlv + muon2_tlv")

    # and compute the invariant mass and momentum of the dimuon system
    df = df.Define("zmumu_m", "zmumu.M()")
    df = df.Define("zmumu_p", "zmumu.P()")


    #########
    ### CUT 2: Z mass window
    #########
    hists.append(df.Histo1D(("zmumu_m", "", *bins_m), "zmumu_m"))
    df = df.Filter("zmumu_m > 86 && zmumu_m < 96")
    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut2"))


    #########
    ### CUT 3: Z momentum
    #########
    hists.append(df.Histo1D(("zmumu_p", "", *bins_m), "zmumu_p"))
    df = df.Filter("zmumu_p > 20 && zmumu_p < 70")
    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut3"))





    #########
    ### CUT 4: cosThetaMiss
    #########
    df = df.Define("missingEnergy", "FCCAnalyses::missingEnergy(240, ReconstructedParticles)")
    df = df.Define("cosTheta_miss", "FCCAnalyses::get_cosTheta_miss(missingEnergy)")
    hists.append(df.Histo1D(("cosThetaMiss", "", *bins_cosThetaMiss), "cosTheta_miss"))
    df = df.Filter("cosTheta_miss < 0.98")
    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut4"))



    #########
    ### CUT 5: recoil cut
    #########

    # compute the recoil
    df = df.Define("init_tlv", "TLorentzVector ret; ret.SetPxPyPzE(0, 0, 0, 240); return ret;") # initial LorentzVector

    df = df.Define("zmumu_recoil_tlv", "init_tlv-zmumu")
    df = df.Define("zmumu_recoil_m", "zmumu_recoil_tlv.M()")


    hists.append(df.Histo1D(("zmumu_recoil_m", "", *bins_m), "zmumu_recoil_m"))
    df = df.Filter("zmumu_recoil_m < 140 && zmumu_recoil_m > 120")
    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut5"))




    # Now we go for jet clustering
    # we want to look for 2 jets in our event, but first we need to remove the muons
    df = df.Define("rps_no_muons", "FCCAnalyses::ReconstructedParticle::remove(ReconstructedParticles, muons)")

    ## define jet clustering parameters
    jetClusteringHelper = ExclusiveJetClusteringHelper("rps_no_muons", 2, "") # force 2 jets clustering
    df = jetClusteringHelper.define(df) # run jet clustering

    ## define jet flavour tagging parameters

    jetFlavourHelper = JetFlavourHelper(
        collections,
        jetClusteringHelper.jets,
        jetClusteringHelper.constituents
    )

    ## define observables for tagger
    df = jetFlavourHelper.define(df)

    ## tagger inference
    df = jetFlavourHelper.inference(weaver_preproc, weaver_model, df)

    ###
    hists.append(df.Histo1D(("recojet_isB", "", *(100, 0, 1)), "recojet_isB"))
    hists.append(df.Histo1D(("recojet_isC", "", *(100, 0, 1)), "recojet_isC"))
    hists.append(df.Histo1D(("recojet_isS", "", *(100, 0, 1)), "recojet_isC"))
    # "recojet_isS", "recojet_isG", "recojet_isU", "recojet_isD", "recojet_isTAU"

    df = df.Define("jet1", f"{jetClusteringHelper.jets}[0]")
    df = df.Define("jet2", f"{jetClusteringHelper.jets}[1]")

    # extract jet kinematics and make LorentzVectors
    df = df.Define("jets_e", f"FCCAnalyses::JetClusteringUtils::get_e({jetClusteringHelper.jets})")
    df = df.Define("jets_px", f"FCCAnalyses::JetClusteringUtils::get_px({jetClusteringHelper.jets})")
    df = df.Define("jets_py", f"FCCAnalyses::JetClusteringUtils::get_py({jetClusteringHelper.jets})")
    df = df.Define("jets_pz", f"FCCAnalyses::JetClusteringUtils::get_pz({jetClusteringHelper.jets})")
    df = df.Define("jet1_tlv", "ROOT::Math::PxPyPzEVector(jets_px[0], jets_py[0], jets_pz[0], jets_e[0])")
    df = df.Define("jet2_tlv", "ROOT::Math::PxPyPzEVector(jets_px[1], jets_py[1], jets_pz[1], jets_e[1])")
    df = df.Define("dijet_tlv", "jet1_tlv + jet2_tlv")
    df = df.Define("dijet_m", "dijet_tlv.M()")


    hists.append(df.Histo1D(("dijet_m", "", *bins_m), "dijet_m"))

    # cut on the b-jet score
    df = df.Filter("recojet_isB[0] > 0.5 && recojet_isB[1] > 0.5")
    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut6"))

    hists.append(df.Histo1D(("zmumu_recoil_m_final", "", *(40, 120, 140)), "zmumu_recoil_m"))

    return hists, weightsum
