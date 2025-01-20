
import ROOT
ROOT.TH1.SetDefaultSumw2(ROOT.kTRUE)



# list of all processes
fraction = 1
processList = {
    'wzp6_ee_ssH_HWW_ecm240':       {'fraction':1},
    'wzp6_ee_ccH_HWW_ecm240':       {'fraction':1},
    'wzp6_ee_bbH_HWW_ecm240':       {'fraction':1},
    'wzp6_ee_qqH_HWW_ecm240':       {'fraction':1},
    'wzp6_ee_nunuH_HWW_ecm240':     {'fraction':1},
    'wzp6_ee_eeH_HWW_ecm240':       {'fraction':1},
    'wzp6_ee_mumuH_HWW_ecm240':     {'fraction':1},
    'wzp6_ee_tautauH_HWW_ecm240':   {'fraction':1},
}


inputDir = "/ceph/submit/data/group/fcc/ee/generation/DelphesEvents/winter2023/IDEA/"
procDict = "/ceph/submit/data/group/fcc/ee/generation/DelphesEvents/winter2023/IDEA/samplesDict.json"

# additional/custom C++ functions
includePaths = ["../../functions/functions.h", "../../functions/functions_gen.h"]


# output directory
outputDir   = "output/h_ww/histmaker/"


# optional: ncpus, default is 4, -1 uses all cores available
nCPUS       = 32

# scale the histograms with the cross-section and integrated luminosity
doScale = True
intLumi = 10800000

# define histograms
bins_m = (250, 0, 250) # 100 MeV bins





def build_graph(df, dataset):

    hists, cols = [], []

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


    # define collections
    df = df.Alias("Particle0", "Particle#0.index")
    df = df.Alias("Particle1", "Particle#1.index")
    df = df.Alias("MCRecoAssociations0", "MCRecoAssociations#0.index")
    df = df.Alias("MCRecoAssociations1", "MCRecoAssociations#1.index")

    # define muons
    df = df.Alias("Muon0", "Muon#0.index")
    df = df.Define("muons_all", "FCCAnalyses::ReconstructedParticle::get(Muon0, ReconstructedParticles)")
    df = df.Define("muons_all_p", "FCCAnalyses::ReconstructedParticle::get_p(muons_all)")
    hists.append(df.Histo1D(("muons_all_p", "", *bins_m), "muons_all_p"))

    df = df.Define("muons_hard", "FCCAnalyses::ReconstructedParticle::sel_p(25)(muons_all)")
    df = df.Define("muons_hard_tlv", "FCCAnalyses::makeLorentzVectors(muons_hard)")
    df = df.Define("muons_p", "FCCAnalyses::ReconstructedParticle::get_p(muons_hard)")
    df = df.Define("muons_no", "FCCAnalyses::ReconstructedParticle::get_n(muons_hard)")
    df = df.Define("muons_q", "FCCAnalyses::ReconstructedParticle::get_charge(muons_hard)")

    df = df.Define("muons_soft", "FCCAnalyses::ReconstructedParticle::sel_p(10)(muons_all)")
    df = df.Define("muons_soft_tlv", "FCCAnalyses::makeLorentzVectors(muons_soft)")
    df = df.Define("muons_soft_p", "FCCAnalyses::ReconstructedParticle::get_p(muons_soft)")
    df = df.Define("muons_soft_no", "FCCAnalyses::ReconstructedParticle::get_n(muons_soft)")
    df = df.Define("muons_soft_q", "FCCAnalyses::ReconstructedParticle::get_charge(muons_soft)")

    # define electrons
    df = df.Alias("Electron0", "Electron#0.index")
    df = df.Define("electrons_all", "FCCAnalyses::ReconstructedParticle::get(Electron0, ReconstructedParticles)")
    df = df.Define("electrons_all_p", "FCCAnalyses::ReconstructedParticle::get_p(electrons_all)")
    hists.append(df.Histo1D(("electrons_all_p", "", *bins_m), "electrons_all_p"))

    df = df.Define("electrons_hard", "FCCAnalyses::ReconstructedParticle::sel_p(25)(electrons_all)")
    df = df.Define("electrons_hard_tlv", "FCCAnalyses::makeLorentzVectors(electrons_hard)")
    df = df.Define("electrons_p", "FCCAnalyses::ReconstructedParticle::get_p(electrons_hard)")
    df = df.Define("electrons_no", "FCCAnalyses::ReconstructedParticle::get_n(electrons_hard)")
    df = df.Define("electrons_q", "FCCAnalyses::ReconstructedParticle::get_charge(electrons_hard)")

    df = df.Define("electrons_soft", "FCCAnalyses::ReconstructedParticle::sel_p(10)(electrons_all)")
    df = df.Define("electrons_soft_tlv", "FCCAnalyses::makeLorentzVectors(electrons_soft)")
    df = df.Define("electrons_soft_p", "FCCAnalyses::ReconstructedParticle::get_p(electrons_soft)")
    df = df.Define("electrons_soft_no", "FCCAnalyses::ReconstructedParticle::get_n(electrons_soft)")
    df = df.Define("electrons_soft_soft_q", "FCCAnalyses::ReconstructedParticle::get_charge(electrons_soft)")

    # missing energy/mass
    df = df.Define("missingMass", "FCCAnalyses::missingMass(240., ReconstructedParticles)")
    hists.append(df.Histo1D(("missingMass", "", *bins_m), "missingMass"))


    return hists, weightsum
