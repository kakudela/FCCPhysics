
import ROOT
ROOT.TH1.SetDefaultSumw2(ROOT.kTRUE)



# list of all processes
fraction = 0.05
processList = {
    'wzp6_ee_mumuH_ecm240':     {'fraction':1, "crossSection": 0.0104},
    'wzp6_ee_mumu_ecm240':     {'fraction':1},
    'wzp6_ee_tautau_ecm240':     {'fraction':1},
}





inputDir = "/ceph/submit/data/group/fcc/ee/generation/DelphesEvents/winter2023/IDEA/"
procDict = "/ceph/submit/data/group/fcc/ee/generation/DelphesEvents/winter2023/IDEA/samplesDict.json"

# additional/custom C++ functions
includePaths = ["../../functions/functions.h", "../../functions/functions_gen.h"]


# output directory
outputDir   = "output/h_zh/histmaker/"


# optional: ncpus, default is 4, -1 uses all cores available
nCPUS       = 32

# scale the histograms with the cross-section and integrated luminosity
doScale = True
intLumi = 250000

# define histograms
bins_m = (250, 0, 250) # 100 MeV bins
bins_maa = (100, 120, 130) # 100 MeV bins
bins_p = (200, 0, 200) # 100 MeV bins
bins_m_zoom = (100, 120, 130) # 100 MeV

bins_theta = (500, 0, 5)
bins_phi = (400, -4, 4)

bins_count = (100, 0, 100)
bins_pdgid = (60, -30, 30)
bins_charge = (10, -5, 5)

bins_resolution = (10000, 0.95, 1.05)
bins_resolution_1 = (20000, 0, 2)

jet_energy = (1000, 0, 100) # 100 MeV bins
dijet_m = (2000, 0, 200) # 100 MeV bins
visMass = (2000, 0, 200) # 100 MeV bins
missEnergy  = (2000, 0, 200) # 100 MeV bins

dijet_m_final = (500, 50, 100) # 100 MeV bins

bins_cos = (100, -1, 1)
bins_cos_abs = (100, 0, 1)
bins_iso = (1000, 0, 10)
bins_aco = (1000,0,5)
bins_cosThetaMiss = (10000, 0, 1)

bins_m_fine = (500, 110, 130) # 100 MeV bins


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

    df = df.Define("ecm", "240")
    df = df.Alias("Lepton0", "Muon#0.index")
    df = df.Define("leps_all", "FCCAnalyses::ReconstructedParticle::get(Lepton0, ReconstructedParticles)")

    # all leptons (bare)
    df = df.Define("leps_all_p", "FCCAnalyses::ReconstructedParticle::get_p(leps_all)")
    df = df.Define("leps_all_theta", "FCCAnalyses::ReconstructedParticle::get_theta(leps_all)")
    df = df.Define("leps_all_phi", "FCCAnalyses::ReconstructedParticle::get_phi(leps_all)")
    df = df.Define("leps_all_q", "FCCAnalyses::ReconstructedParticle::get_charge(leps_all)")
    df = df.Define("leps_all_no", "FCCAnalyses::ReconstructedParticle::get_n(leps_all)")
    df = df.Define("leps_all_iso", "FCCAnalyses::coneIsolation(0.01, 0.5)(leps_all, ReconstructedParticles)") 

    # cuts on leptons
    df = df.Define("leps", "FCCAnalyses::ReconstructedParticle::sel_p(20)(leps_all)")

    df = df.Define("leps_p", "FCCAnalyses::ReconstructedParticle::get_p(leps)")
    df = df.Define("leps_theta", "FCCAnalyses::ReconstructedParticle::get_theta(leps)")
    df = df.Define("leps_phi", "FCCAnalyses::ReconstructedParticle::get_phi(leps)")
    df = df.Define("leps_q", "FCCAnalyses::ReconstructedParticle::get_charge(leps)")
    df = df.Define("leps_no", "FCCAnalyses::ReconstructedParticle::get_n(leps)")
    df = df.Define("leps_iso", "FCCAnalyses::coneIsolation(0.01, 0.5)(leps, ReconstructedParticles)")
    #df = df.Define("leps_sel_iso", "FCCAnalyses::sel_iso(0.25)(leps, leps_iso)") # 0.25
    df = df.Define("leps_sel_iso", "FCCAnalyses::sel_range(0, 0.25, false)(leps, leps_iso)") # based on vvH sample


    # baseline selections and histograms
    hists.append(df.Histo1D(("leps_all_p_cut0", "", *bins_m), "leps_all_p"))
    hists.append(df.Histo1D(("leps_all_theta_cut0", "", *bins_theta), "leps_all_theta"))
    hists.append(df.Histo1D(("leps_all_phi_cut0", "", *bins_phi), "leps_all_phi"))
    hists.append(df.Histo1D(("leps_all_q_cut0", "", *bins_charge), "leps_all_q"))
    hists.append(df.Histo1D(("leps_all_no_cut0", "", *bins_count), "leps_all_no"))
    hists.append(df.Histo1D(("leps_all_iso_cut0", "", *bins_iso), "leps_all_iso"))


    hists.append(df.Histo1D(("leps_p_cut0", "", *bins_m), "leps_p"))
    hists.append(df.Histo1D(("leps_theta_cut0", "", *bins_theta), "leps_theta"))
    hists.append(df.Histo1D(("leps_phi_cut0", "", *bins_phi), "leps_phi"))
    hists.append(df.Histo1D(("leps_q_cut0", "", *bins_charge), "leps_q"))
    hists.append(df.Histo1D(("leps_no_cut0", "", *bins_count), "leps_no"))
    hists.append(df.Histo1D(("leps_iso_cut0", "", *bins_iso), "leps_iso"))



    
    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut0"))

    #########
    ### CUT 1: at least a lepton with at least 1 isolated one
    #########
    df = df.Filter("leps_no >= 1 && leps_sel_iso.size() > 0")
    hists.append(df.Histo1D(("cutFlow_cut1", "", *bins_count), "cut1"))
    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut1"))




    #########
    ### CUT 2 :at least 2 OS leptons, and build the resonance
    #########
    df = df.Filter("leps_no >= 2 && abs(Sum(leps_q)) < leps_q.size()")
    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut2"))

    #df = df.Filter("leps_no == 2")

    # build the Z resonance based on the available leptons. Returns the best lepton pair compatible with the Z mass and recoil at 125 GeV
    # technically, it returns a ReconstructedParticleData object with index 0 the di-lepton system, index and 2 the leptons of the pair
    df = df.Define("zbuilder_result", "FCCAnalyses::resonanceBuilder_mass_recoil(91.2, 125, 0.4, ecm, false)(leps, MCRecoAssociations0, MCRecoAssociations1, ReconstructedParticles, Particle, Particle0, Particle1)")
    df = df.Define("zll", "ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>{zbuilder_result[0]}") # the Z
    df = df.Define("zll_leps", "ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>{zbuilder_result[1],zbuilder_result[2]}") # the leptons
    df = df.Define("zll_m", "FCCAnalyses::ReconstructedParticle::get_mass(zll)[0]")
    df = df.Define("zll_p", "FCCAnalyses::ReconstructedParticle::get_p(zll)[0]")
    df = df.Define("zll_recoil", "FCCAnalyses::ReconstructedParticle::recoilBuilder(ecm)(zll)")
    df = df.Define("zll_recoil_m", "FCCAnalyses::ReconstructedParticle::get_mass(zll_recoil)[0]")
    df = df.Define("zll_leps_p", "FCCAnalyses::ReconstructedParticle::get_p(zll_leps)")
    df = df.Define("zll_leps_theta", "FCCAnalyses::ReconstructedParticle::get_theta(zll_leps)")


    df = df.Define("missingEnergy", "FCCAnalyses::missingEnergy(ecm, ReconstructedParticles)")
    df = df.Define("cosTheta_miss", "FCCAnalyses::get_cosTheta_miss(missingEnergy)")
    #df = df.Define("cosTheta_miss", "FCCAnalyses::get_cosTheta_miss(MissingET)")

    df = df.Define("acoplanarity", "FCCAnalyses::acoplanarity(leps)")
    df = df.Define("acolinearity", "FCCAnalyses::acolinearity(leps)")
    
    hists.append(df.Histo1D(("zll_m_cut2", "", *bins_m), "zll_m"))


    #########
    ### CUT 3: Z mass window
    #########  
    df = df.Filter("zll_m > 86 && zll_m < 96")
    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut3"))
    hists.append(df.Histo1D(("zll_p_cut3", "", *bins_m), "zll_p"))

    #########
    ### CUT 4: Z momentum
    #########  
    df = df.Filter("zll_p > 20 && zll_p < 70")
    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut4"))
    hists.append(df.Histo1D(("zll_recoil_cut4", "", *bins_m), "zll_recoil_m"))

    #########
    ### CUT 5: recoil cut
    #########  
    df = df.Filter("zll_recoil_m > 110 && zll_recoil_m < 150")
    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut5"))
    hists.append(df.Histo1D(("cosThetaMiss_cut5", "", *bins_cosThetaMiss), "cosTheta_miss"))

    #########
    ### CUT 6: cosThetaMiss, for mass analysis
    #########  
    df = df.Filter("cosTheta_miss < 0.98")

    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut6"))

    # final histograms
    hists.append(df.Histo1D(("leps_p", "", *bins_m), "leps_p"))
    hists.append(df.Histo1D(("zll_p", "", *bins_m), "zll_p"))
    hists.append(df.Histo1D(("zll_m", "", *bins_m), "zll_m"))
    hists.append(df.Histo1D(("zll_recoil", "", *(40, 110, 150)), "zll_recoil_m"))

    hists.append(df.Histo1D(("cosThetaMiss", "", *bins_cosThetaMiss), "cosTheta_miss"))
    hists.append(df.Histo1D(("acoplanarity", "", *bins_aco), "acoplanarity"))
    hists.append(df.Histo1D(("acolinearity", "", *bins_aco), "acolinearity"))


    return hists, weightsum
