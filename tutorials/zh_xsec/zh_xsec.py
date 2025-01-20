
import ROOT
ROOT.TH1.SetDefaultSumw2(ROOT.kTRUE)



# list of all processes
processList = {
    'wzp6_ee_mumuH_ecm240': {'fraction':1},
    'p8_ee_ZZ_ecm240': {'fraction':1},
    'p8_ee_WW_mumu_ecm240': {'fraction':1},
    'wzp6_ee_mumu_ecm240': {'fraction':1},
}



inputDir = "/ceph/submit/data/group/fcc/ee/generation/DelphesEvents/winter2023/IDEA/"
procDict = "/ceph/submit/data/group/fcc/ee/generation/DelphesEvents/winter2023/IDEA/samplesDict.json"

# additional/custom C++ functions
includePaths = ["../../functions/functions.h", "../../functions/functions_gen.h"]


# output directory
outputDir   = "output/tutorials/zh_xsec/"


# optional: ncpus, default is 4, -1 uses all cores available
nCPUS       = 32

# scale the histograms with the cross-section and integrated luminosity
doScale = True
intLumi = 10800000 # 10.8 ab-1

# define histograms
bins_p_mu = (250, 0, 250)
bins_m_ll = (250, 0, 250)
bins_p_ll = (250, 0, 250)
bins_recoil = (250, 0, 250)
bins_cosThetaMiss = (10000, 0, 1)
bins_recoil_final = (200, 120, 140)

bins_theta = (500, -5, 5)
bins_eta = (600, -3, 3)
bins_phi = (500, -5, 5)

bins_count = (50, 0, 50)
bins_charge = (10, -5, 5)
bins_iso = (500, 0, 5)



def build_graph(df, dataset):

    hists = []
    df = df.Define("weight", "1.0")
    weightsum = df.Sum("weight")

    # define some aliases to be used later on
    df = df.Alias("Particle0", "Particle#0.index")
    df = df.Alias("Particle1", "Particle#1.index")
    df = df.Alias("MCRecoAssociations0", "MCRecoAssociations#0.index")
    df = df.Alias("MCRecoAssociations1", "MCRecoAssociations#1.index")
    df = df.Alias("Muon0", "Muon#0.index")


    # get all the leptons from the collection
    df = df.Define("muons_all", "FCCAnalyses::ReconstructedParticle::get(Muon0, ReconstructedParticles)")


    # select leptons with momentum > 20 GeV
    df = df.Define("muons", "FCCAnalyses::ReconstructedParticle::sel_p(20)(muons_all)")

    df = df.Define("muons_p", "FCCAnalyses::ReconstructedParticle::get_p(muons)")
    df = df.Define("muons_theta", "FCCAnalyses::ReconstructedParticle::get_theta(muons)")
    df = df.Define("muons_phi", "FCCAnalyses::ReconstructedParticle::get_phi(muons)")
    df = df.Define("muons_q", "FCCAnalyses::ReconstructedParticle::get_charge(muons)")
    df = df.Define("muons_no", "FCCAnalyses::ReconstructedParticle::get_n(muons)")

    # compute the muon isolation and store muons with an isolation cut of 0.25 in a separate column muons_sel_iso
    df = df.Define("muons_iso", "FCCAnalyses::coneIsolation(0.01, 0.5)(muons, ReconstructedParticles)")
    df = df.Define("muons_sel_iso", "FCCAnalyses::sel_range(0, 0.25, false)(muons, muons_iso)")

    # baseline histograms, before any selection cuts (store with _cut0)
    hists.append(df.Histo1D(("muons_p_cut0", "", *bins_p_mu), "muons_p"))
    hists.append(df.Histo1D(("muons_theta_cut0", "", *bins_theta), "muons_theta"))
    hists.append(df.Histo1D(("muons_phi_cut0", "", *bins_phi), "muons_phi"))
    hists.append(df.Histo1D(("muons_q_cut0", "", *bins_charge), "muons_q"))
    hists.append(df.Histo1D(("muons_no_cut0", "", *bins_count), "muons_no"))
    hists.append(df.Histo1D(("muons_iso_cut0", "", *bins_iso), "muons_iso"))

    #########
    ### CUT 0: all events
    #########
    df = df.Define("cut0", "0")
    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut0"))


    #########
    ### CUT 1: at least 1 muon with at least one isolated one
    #########
    df = df.Filter("muons_no >= 1 && muons_sel_iso.size() > 0")
    df = df.Define("cut1", "1")
    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut1"))


    #########
    ### CUT 2 :at least 2 opposite-sign (OS) leptons
    #########
    df = df.Filter("muons_no >= 2 && abs(Sum(muons_q)) < muons_q.size()")
    df = df.Define("cut2", "2")
    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut2"))

    # now we build the Z resonance based on the available leptons.
    # the function resonanceBuilder_mass_recoil returns the best lepton pair compatible with the Z mass (91.2 GeV) and recoil at 125 GeV
    # the argument 0.4 gives a weight to the Z mass and the recoil mass in the chi2 minimization
    # technically, it returns a ReconstructedParticleData object with index 0 the di-lepton system, index and 2 the leptons of the pair
    df = df.Define("zbuilder_result", "FCCAnalyses::resonanceBuilder_mass_recoil(91.2, 125, 0.4, 240, false)(muons, MCRecoAssociations0, MCRecoAssociations1, ReconstructedParticles, Particle, Particle0, Particle1)")
    df = df.Define("zmumu", "Vec_rp{zbuilder_result[0]}") # the Z
    df = df.Define("zmumu_muons", "Vec_rp{zbuilder_result[1],zbuilder_result[2]}") # the leptons 
    df = df.Define("zmumu_m", "FCCAnalyses::ReconstructedParticle::get_mass(zmumu)[0]") # Z mass
    df = df.Define("zmumu_p", "FCCAnalyses::ReconstructedParticle::get_p(zmumu)[0]") # momentum of the Z
    df = df.Define("zmumu_recoil", "FCCAnalyses::ReconstructedParticle::recoilBuilder(240)(zmumu)") # compute the recoil based on the reconstructed Z
    df = df.Define("zmumu_recoil_m", "FCCAnalyses::ReconstructedParticle::get_mass(zmumu_recoil)[0]") # recoil mass
    df = df.Define("zmumu_muons_p", "FCCAnalyses::ReconstructedParticle::get_p(zmumu_muons)") # get the momentum of the 2 muons from the Z resonance


    #########
    ### CUT 3: Z mass window
    #########
    hists.append(df.Histo1D(("zmumu_m", "", *bins_m_ll), "zmumu_m"))
    df = df.Filter("zmumu_m > 86 && zmumu_m < 96")
    df = df.Define("cut3", "3")
    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut3"))

    #########
    ### CUT 4: Z momentum
    #########
    hists.append(df.Histo1D(("zmumu_p", "", *bins_p_ll), "zmumu_p"))
    df = df.Filter("zmumu_p > 20 && zmumu_p < 70")
    df = df.Define("cut4", "4")
    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut4"))


    #########
    ### CUT 5: cosThetaMiss
    #########
    df = df.Define("missingEnergy", "FCCAnalyses::missingEnergy(240., ReconstructedParticles)")
    df = df.Define("cosTheta_miss", "FCCAnalyses::get_cosTheta_miss(missingEnergy)")
    hists.append(df.Histo1D(("cosThetaMiss_cut4", "", *bins_cosThetaMiss), "cosTheta_miss")) # plot it before the cut
    df = df.Filter("cosTheta_miss < 0.98")
    df = df.Define("cut5", "5")
    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut5"))


    #########
    ### CUT 6: recoil mass window
    #########
    hists.append(df.Histo1D(("zmumu_recoil_m", "", *bins_recoil), "zmumu_recoil_m"))
    df = df.Filter("zmumu_recoil_m < 140 && zmumu_recoil_m > 120")
    df = df.Define("cut6", "6")
    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut6"))
    hists.append(df.Histo1D(("zmumu_recoil_m_final", "", *bins_recoil_final), "zmumu_recoil_m"))


    return hists, weightsum
