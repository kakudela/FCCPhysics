
import ROOT
ROOT.TH1.SetDefaultSumw2(ROOT.kTRUE)



# list of all processes
fraction = 1
processList = {
    'wzp6_ee_nunuH_Haa_ecm240':     {'fraction':1},
    'wzp6_ee_eeH_Haa_ecm240':       {'fraction':1},
    'wzp6_ee_tautauH_Haa_ecm240':   {'fraction':1},
    'wzp6_ee_ccH_Haa_ecm240':       {'fraction':1},
    'wzp6_ee_bbH_Haa_ecm240':       {'fraction':1},
    'wzp6_ee_qqH_Haa_ecm240':       {'fraction':1},
    'wzp6_ee_ssH_Haa_ecm240':       {'fraction':1},
    'wzp6_ee_mumuH_Haa_ecm240':     {'fraction':1},
    #'kkmcee_ee_uu_ecm240':          {'fraction':fraction},
    #'kkmcee_ee_dd_ecm240':          {'fraction':fraction},
    #'kkmcee_ee_cc_ecm240':          {'fraction':fraction},
    #'kkmcee_ee_ss_ecm240':          {'fraction':fraction},
    #'kkmcee_ee_bb_ecm240':          {'fraction':fraction},
    #'kkmcee_ee_tautau_ecm240':      {'fraction':fraction},
    #'kkmcee_ee_mumu_ecm240':        {'fraction':fraction},
    #'kkmcee_ee_nuenue_ecm240':      {'fraction':fraction},
    #'kkmcee_ee_numunumu_ecm240':    {'fraction':fraction},
    #'kkmcee_ee_nutaunutau_ecm240':  {'fraction':fraction},
    'wzp6_ee_gammagamma_ecm240':    {'fraction':fraction},
    #'yfsww_ee_ww_ecm240':           {'fraction':fraction},
    'wz3p6_ee_uu_ecm240':           {'fraction':fraction},
    'wz3p6_ee_dd_ecm240':           {'fraction':fraction},
    'wz3p6_ee_cc_ecm240':           {'fraction':fraction},
    'wz3p6_ee_ss_ecm240':           {'fraction':fraction},
    'wz3p6_ee_bb_ecm240':           {'fraction':fraction},
    'wz3p6_ee_tautau_ecm240':       {'fraction':fraction},
    'wz3p6_ee_mumu_ecm240':         {'fraction':fraction},
    'wz3p6_ee_ee_Mee_30p_ecm240':   {'fraction':fraction},
    'wz3p6_ee_nunu_ecm240':         {'fraction':fraction},
}





inputDir = "/ceph/submit/data/group/fcc/ee/generation/DelphesEvents/winter2023/IDEA/"
procDict = "/ceph/submit/data/group/fcc/ee/generation/DelphesEvents/winter2023/IDEA/samplesDict.json"

# additional/custom C++ functions
includePaths = ["../../functions/functions.h", "../../functions/functions_gen.h"]


# output directory
outputDir   = "output/h_aa/histmaker/"


# optional: ncpus, default is 4, -1 uses all cores available
nCPUS       = 32

# scale the histograms with the cross-section and integrated luminosity
doScale = True
intLumi = 10800000

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

    # all photons
    df = df.Alias("Photon0", "Photon#0.index")
    df = df.Define("photons_all", "FCCAnalyses::ReconstructedParticle::get(Photon0, ReconstructedParticles)")
    df = df.Define("photons_all_p", "FCCAnalyses::ReconstructedParticle::get_p(photons_all)")
    df = df.Define("photons_all_theta", "FCCAnalyses::ReconstructedParticle::get_theta(photons_all)")
    df = df.Define("photons_all_no", "FCCAnalyses::ReconstructedParticle::get_n(photons_all)")
    df = df.Define("photons_all_costheta", "FCCAnalyses::Vec_f ret; for(auto & theta: photons_all_theta) ret.push_back(std::abs(cos(theta))); return ret;")

    hists.append(df.Histo1D(("photons_all_p", "", *bins_p), "photons_all_p"))
    hists.append(df.Histo1D(("photons_all_costheta", "", *bins_cos_abs), "photons_all_costheta"))

    # select photon momentum 40-95 GeV
    df = df.Define("photons", "FCCAnalyses::sel_range(40, 95, false)(photons_all, photons_all_p)")
    df = df.Define("photons_p", "FCCAnalyses::ReconstructedParticle::get_p(photons)")
    df = df.Define("photons_theta", "FCCAnalyses::ReconstructedParticle::get_theta(photons)")
    df = df.Define("photons_n", "FCCAnalyses::ReconstructedParticle::get_n(photons)")
    df = df.Define("photons_costheta", "FCCAnalyses::Vec_f ret; for(auto & theta: photons_theta) ret.push_back(std::abs(cos(theta))); return ret;")
    df = df.Define("photons_iso", "FCCAnalyses::coneIsolation(0.01, 0.5)(photons, ReconstructedParticles)")


    #########
    ### CUT 0: all events
    #########
    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut0"))
    hists.append(df.Histo1D(("photons_n", "", *bins_count), "photons_n"))
    hists.append(df.Histo1D(("photons_p", "", *bins_p), "photons_p"))


    #########
    ### CUT 1: at least 1 photon
    #########
    hists.append(df.Histo1D(("photons_n_nOne", "", *bins_count), "photons_n"))
    df = df.Filter("photons_n >= 1") ### ORIG

    df = df.Define("photon_leading_p", "photons_p[0]")
    df = df.Define("photon_leading_costheta", "photons_costheta[0]")
    df = df.Define("photon_leading_iso", "photons_iso[0]")
    hists.append(df.Histo1D(("photon_leading_p", "", *bins_p), "photon_leading_p"))
    hists.append(df.Histo1D(("photon_leading_costheta", "", *bins_cos_abs), "photon_leading_costheta"))
    hists.append(df.Histo1D(("photon_leading_iso", "", *bins_iso), "photon_leading_iso"))

    ####df = df.Filter("photon_leading_iso < 0.5 && photon_leading_costheta < 0.9") ### ORIG
    df = df.Filter("photon_leading_iso < 0.1 && photon_leading_costheta < 0.75") ## new
    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut1"))

    #########
    ### CUT 2: at least 2 photons
    #########
    df = df.Filter("photons_n >= 2") ### ORIG

    df = df.Define("photon_subleading_p", "photons_p[1]")
    df = df.Define("photon_subleading_costheta", "photons_costheta[1]")
    df = df.Define("photon_subleading_iso", "photons_iso[1]")
    hists.append(df.Histo1D(("photon_subleading_p", "", *bins_p), "photon_subleading_p"))
    hists.append(df.Histo1D(("photon_subleading_costheta", "", *bins_cos_abs), "photon_subleading_costheta"))
    hists.append(df.Histo1D(("photon_subleading_iso", "", *bins_iso), "photon_subleading_iso"))

    df = df.Filter("photon_subleading_costheta < 0.85") ## new

    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut2"))

    #########
    ### CUT 3 : exactly 2 photons
    #########
    df = df.Filter("photons_n == 2") ### ORIG
    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut3"))
    df = df.Alias("haa_photons", "photons")

    df = df.Define("photon1_tlv", "ROOT::Math::PxPyPzEVector(haa_photons[0].momentum.x, haa_photons[0].momentum.y, haa_photons[0].momentum.z, haa_photons[0].energy)")
    df = df.Define("photon2_tlv", "ROOT::Math::PxPyPzEVector(haa_photons[1].momentum.x, haa_photons[1].momentum.y, haa_photons[1].momentum.z, haa_photons[1].energy)")
    df = df.Define("haa_tlv", "photon1_tlv+photon2_tlv")
    df = df.Define("haa_m", "haa_tlv.M()")
    df = df.Define("haa_p", "haa_tlv.P()")
    df = df.Define("init_tlv", "ROOT::Math::PxPyPzEVector(0, 0, 0, 240)")
    df = df.Define("haa_recoil_tlv", "init_tlv-haa_tlv")
    df = df.Define("haa_recoil_m", "haa_recoil_tlv.M()")



    #########
    ### CUT 4: recoil cut (Z mass)
    #########  
    hists.append(df.Histo1D(("haa_recoil_m_nOne", "", *bins_m), "haa_recoil_m"))
    df = df.Filter("haa_recoil_m > 85 && haa_recoil_m < 110") # see significance test
    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut4"))


    #####
    ### CUT 5: momentum
    #####
    hists.append(df.Histo1D(("haa_p_nOne", "", *bins_p), "haa_p"))
    #df = df.Filter("haa_p > 30 && haa_p < 55") ## orig
    df = df.Filter("haa_p > 20 && haa_p < 55")
    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut5"))


    df = df.Define("acoplanarity", "FCCAnalyses::acoplanarity(haa_photons)")
    df = df.Define("acolinearity", "FCCAnalyses::acolinearity(haa_photons)")

    ####
    ## CUT 6: acolinearity
    ####
    hists.append(df.Histo1D(("acolinearity_nOne", "", *bins_aco), "acolinearity"))
    ##df = df.Filter("acolinearity > 0.1 && acolinearity < 0.8") # orig
    df = df.Filter("acolinearity > 0.2 && acolinearity < 0.8")
    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut6"))


    ####
    ## CUT 7: acoplanarity
    ####
    hists.append(df.Histo1D(("acoplanarity_nOne", "", *bins_aco), "acoplanarity"))
    ##df = df.Filter("acoplanarity > 0.02") # orig
    df = df.Filter("acoplanarity > 0.05 && acoplanarity < 1.2")
    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut7"))


    ####
    ## CUT 8: cos theta(miss)
    ####
    df = df.Define("missingEnergy_rp", "FCCAnalyses::missingEnergy(240., ReconstructedParticles)")
    df = df.Define("missingEnergy", "missingEnergy_rp[0].energy")
    df = df.Define("cosThetaMiss", "FCCAnalyses::get_cosTheta_miss(missingEnergy_rp)")
    hists.append(df.Histo1D(("cosThetaMiss_nOne", "", *bins_cosThetaMiss), "cosThetaMiss"))
    df = df.Filter("cosThetaMiss < .995") # 0.98
    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut8"))


    #########
    ### CUT 9 :cut on Higgs mass
    #########
    hists.append(df.Histo1D(("haa_m_nOne", "", *bins_maa), "haa_m"))
    df = df.Filter("haa_m > 120 && haa_m < 130")
    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut9"))


    #########
    ### extra variables
    #########
    
    # missing energy
    hists.append(df.Histo1D(("missingEnergy", "", *bins_m), "missingEnergy"))
    
    # missing mass
    df = df.Define("missingMass", "FCCAnalyses::missingMass(240., ReconstructedParticles)")
    hists.append(df.Histo1D(("missingMass", "", *bins_m), "missingMass"))

    # cos theta of diphoton system
    df = df.Define("haa_theta", "haa_tlv.Theta()")
    df = df.Define("haa_costheta", "abs(cos(haa_theta))")
    hists.append(df.Histo1D(("haa_costheta", "", *(100, 0, 1)), "haa_costheta"))
    
    # opening angle between photons
    df = df.Define("vec3_photon1", "ROOT::Math::XYZVector(photon1_tlv.Px(),photon1_tlv.Py(),photon1_tlv.Pz());")
    df = df.Define("vec3_photon2", "ROOT::Math::XYZVector(photon2_tlv.Px(),photon2_tlv.Py(),photon2_tlv.Pz());")
    df = df.Define("vec3_diphoton", "ROOT::Math::XYZVector(haa_tlv.Px(),haa_tlv.Py(),haa_tlv.Pz());")
    df = df.Define("cosThetaPhotons", "abs(cos((vec3_photon1 - vec3_photon2).Theta()))")
    hists.append(df.Histo1D(("cosThetaPhotons", "", *(100, 0, 1)), "cosThetaPhotons"))

    # photon energy ratio ratio
    df = df.Define("photon_momentum_ratio", "photon_leading_p/photon_subleading_p")
    hists.append(df.Histo1D(("photon_momentum_ratio", "", *(100, 0, 10)), "photon_momentum_ratio"))

    # dphi between photons
    df = df.Define("photon_dphi", "ROOT::Math::VectorUtil::DeltaPhi(photon1_tlv, photon2_tlv)")
    hists.append(df.Histo1D(("photon_dphi", "", *(200, -10, 10)), "photon_dphi"))




    ############################################################################


    # muons
    df = df.Alias("Muon0", "Muon#0.index")
    df = df.Define("muons_all", "FCCAnalyses::ReconstructedParticle::get(Muon0, ReconstructedParticles)")
    df = df.Define("muons_all_p", "FCCAnalyses::ReconstructedParticle::get_p(muons_all)")
    df = df.Define("muons_all_theta", "FCCAnalyses::ReconstructedParticle::get_theta(muons_all)")
    df = df.Define("muons_all_phi", "FCCAnalyses::ReconstructedParticle::get_phi(muons_all)")
    df = df.Define("muons_all_q", "FCCAnalyses::ReconstructedParticle::get_charge(muons_all)")
    df = df.Define("muons_all_no", "FCCAnalyses::ReconstructedParticle::get_n(muons_all)")

    df = df.Define("muons", "FCCAnalyses::ReconstructedParticle::sel_p(25)(muons_all)")
    df = df.Define("muons_p", "FCCAnalyses::ReconstructedParticle::get_p(muons)")
    df = df.Define("muons_theta", "FCCAnalyses::ReconstructedParticle::get_theta(muons)")
    df = df.Define("muons_phi", "FCCAnalyses::ReconstructedParticle::get_phi(muons)")
    df = df.Define("muons_q", "FCCAnalyses::ReconstructedParticle::get_charge(muons)")
    df = df.Define("muons_no", "FCCAnalyses::ReconstructedParticle::get_n(muons)")

    
    # electrons
    df = df.Alias("Electron0", "Electron#0.index")
    df = df.Define("electrons_all", "FCCAnalyses::ReconstructedParticle::get(Electron0, ReconstructedParticles)")
    df = df.Define("electrons_all_p", "FCCAnalyses::ReconstructedParticle::get_p(electrons_all)")
    df = df.Define("electrons_all_theta", "FCCAnalyses::ReconstructedParticle::get_theta(electrons_all)")
    df = df.Define("electrons_all_phi", "FCCAnalyses::ReconstructedParticle::get_phi(electrons_all)")
    df = df.Define("electrons_all_q", "FCCAnalyses::ReconstructedParticle::get_charge(electrons_all)")
    df = df.Define("electrons_all_no", "FCCAnalyses::ReconstructedParticle::get_n(electrons_all)")

    df = df.Define("electrons", "FCCAnalyses::ReconstructedParticle::sel_p(25)(electrons_all)")
    df = df.Define("electrons_p", "FCCAnalyses::ReconstructedParticle::get_p(electrons)")
    df = df.Define("electrons_theta", "FCCAnalyses::ReconstructedParticle::get_theta(electrons)")
    df = df.Define("electrons_phi", "FCCAnalyses::ReconstructedParticle::get_phi(electrons)")
    df = df.Define("electrons_q", "FCCAnalyses::ReconstructedParticle::get_charge(electrons)")
    df = df.Define("electrons_no", "FCCAnalyses::ReconstructedParticle::get_n(electrons)")


    # lepton kinematic histograms
    hists.append(df.Histo1D(("muons_all_p", "", *bins_p), "muons_all_p"))
    hists.append(df.Histo1D(("muons_all_theta", "", *bins_theta), "muons_all_theta"))
    hists.append(df.Histo1D(("muons_all_phi", "", *bins_phi), "muons_all_phi"))
    hists.append(df.Histo1D(("muons_all_q", "", *bins_charge), "muons_all_q"))
    hists.append(df.Histo1D(("muons_all_no", "", *bins_count), "muons_all_no"))

    hists.append(df.Histo1D(("electrons_all_p", "", *bins_p), "electrons_all_p"))
    hists.append(df.Histo1D(("electrons_all_theta", "", *bins_theta), "electrons_all_theta"))
    hists.append(df.Histo1D(("electrons_all_phi", "", *bins_phi), "electrons_all_phi"))
    hists.append(df.Histo1D(("electrons_all_q", "", *bins_charge), "electrons_all_q"))
    hists.append(df.Histo1D(("electrons_all_no", "", *bins_count), "electrons_all_no"))



    ##### CATEGORIZATION: based on #muons, # electrons, missing energy
    select_mumu = "muons_no == 2 && electrons_no == 0 && missingMass < 15"
    select_ee = "electrons_no == 2 && muons_no == 0 && missingMass < 15"
    select_nunu = "electrons_no == 0 && muons_no == 0 && missingMass > 85"
    select_qq = "electrons_no == 0 && muons_no == 0 && missingMass < 15"
    #select_tauhtauh = "electrons_no == 0 && muons_no == 0 && missingMass > 15 && missingMass < 85"




    #######
    # qq final state
    #######
    df_qq = df.Filter(select_qq)
    hists.append(df_qq.Histo1D(("cutFlow", "", *bins_count), "cut10"))


    # define PF candidates collection by removing the muons
    df_qq = df_qq.Define("rps_no_muons", "FCCAnalyses::ReconstructedParticle::remove(ReconstructedParticles, photons)")
    df_qq = df_qq.Define("RP_px", "FCCAnalyses::ReconstructedParticle::get_px(rps_no_muons)")
    df_qq = df_qq.Define("RP_py", "FCCAnalyses::ReconstructedParticle::get_py(rps_no_muons)")
    df_qq = df_qq.Define("RP_pz","FCCAnalyses::ReconstructedParticle::get_pz(rps_no_muons)")
    df_qq = df_qq.Define("RP_e", "FCCAnalyses::ReconstructedParticle::get_e(rps_no_muons)")
    df_qq = df_qq.Define("RP_m", "FCCAnalyses::ReconstructedParticle::get_mass(rps_no_muons)")
    df_qq = df_qq.Define("RP_q", "FCCAnalyses::ReconstructedParticle::get_charge(rps_no_muons)")
    df_qq = df_qq.Define("pseudo_jets", "FCCAnalyses::JetClusteringUtils::set_pseudoJets(RP_px, RP_py, RP_pz, RP_e)")

    df_qq = df_qq.Define("clustered_jets", "JetClustering::clustering_ee_kt(2, 2, 1, 0)(pseudo_jets)")
    df_qq = df_qq.Define("jets", "FCCAnalyses::JetClusteringUtils::get_pseudoJets(clustered_jets)")
    df_qq = df_qq.Define("jetconstituents", "FCCAnalyses::JetClusteringUtils::get_constituents(clustered_jets)")
    df_qq = df_qq.Define("jets_e", "FCCAnalyses::JetClusteringUtils::get_e(jets)")
    df_qq = df_qq.Define("jets_px", "FCCAnalyses::JetClusteringUtils::get_px(jets)")
    df_qq = df_qq.Define("jets_py", "FCCAnalyses::JetClusteringUtils::get_py(jets)")
    df_qq = df_qq.Define("jets_pz", "FCCAnalyses::JetClusteringUtils::get_pz(jets)")
    df_qq = df_qq.Define("jets_m", "FCCAnalyses::JetClusteringUtils::get_m(jets)")

    df_qq = df_qq.Define("jet1", "ROOT::Math::PxPyPzEVector(jets_px[0], jets_py[0], jets_pz[0], jets_e[0])")
    df_qq = df_qq.Define("jet2", "ROOT::Math::PxPyPzEVector(jets_px[1], jets_py[1], jets_pz[1], jets_e[1])")
    df_qq = df_qq.Define("jet1_p", "jet1.P()")
    df_qq = df_qq.Define("jet2_p", "jet2.P()")
    df_qq = df_qq.Define("jet1_theta", "jet1.Theta()")
    df_qq = df_qq.Define("jet2_theta", "jet2.Theta()")
    df_qq = df_qq.Define("dijet", "jet1+jet2")
    df_qq = df_qq.Define("dijet_tlv", "TLorentzVector ret; ret.SetPxPyPzE(dijet.Px(), dijet.Py(), dijet.Pz(), dijet.E()); return ret;")
    df_qq = df_qq.Define("dijet_m", "dijet.M()")
    df_qq = df_qq.Define("dijet_p", "dijet.P()")
    df_qq = df_qq.Define("dijet_theta", "dijet.Theta()")
    df_qq = df_qq.Define("costheta1", "abs(cos(dijet.Theta()))")

    hists.append(df_qq.Histo1D(("zqq_qq_p_nOne", "", *bins_p), "dijet_p"))
    df_qq = df_qq.Filter("dijet_p > 25 && dijet_p < 55")

    hists.append(df_qq.Histo1D(("zqq_m_nOne", "", *bins_m), "dijet_m"))
    df_qq = df_qq.Filter("dijet_m < 105 && dijet_m > 75")

    hists.append(df_qq.Histo1D(("zqq_m", "", *bins_m), "dijet_m"))
    hists.append(df_qq.Histo1D(("zqq_haa_m", "", *bins_m_zoom), "haa_m"))
    hists.append(df_qq.Histo1D(("zqq_aa_p", "", *bins_p), "haa_p"))
    hists.append(df_qq.Histo1D(("zqq_qq_p", "", *bins_p), "dijet_p"))

    hists.append(df_qq.Histo1D(("zqq_acoplanarity", "", *bins_aco), "acoplanarity"))
    hists.append(df_qq.Histo1D(("zqq_acolinearity", "", *bins_aco), "acolinearity"))
    hists.append(df_qq.Histo1D(("zqq_cosThetaMiss", "", *bins_cosThetaMiss), "cosThetaMiss"))



    # extra variables

    # cos(theta1,2) CP analysis inspired
    df_qq = df_qq.Define("vec3_q1", "ROOT::Math::XYZVector(jet1.Px(),jet1.Py(),jet1.Pz());")
    df_qq = df_qq.Define("vec3_q2", "ROOT::Math::XYZVector(jet2.Px(),jet2.Py(),jet2.Pz());")
    df_qq = df_qq.Define("vec3_qq", "ROOT::Math::XYZVector(dijet.Px(),dijet.Py(),dijet.Pz());")
    df_qq = df_qq.Define("costheta2_1", "abs(cos((vec3_q1-vec3_qq).Theta()))")
    df_qq = df_qq.Define("costheta2_2", "abs(cos((vec3_q2-vec3_qq).Theta()))")

    hists.append(df_qq.Histo1D(("zqq_costheta1", "", *(100, 0, 1)), "costheta1"))
    hists.append(df_qq.Histo1D(("zqq_costheta2_1", "", *(100, 0, 1)), "costheta2_1"))
    hists.append(df_qq.Histo1D(("zqq_costheta2_2", "", *(100, 0, 1)), "costheta2_2"))


    # angle between Higgs and Z
    df_qq = df_qq.Define("photon_qq_costheta", "abs(cos((vec3_diphoton-vec3_qq).Theta()))")
    hists.append(df_qq.Histo1D(("zqq_photon_qq_costheta", "", *(100, -1, 1)), "photon_qq_costheta"))


    #######
    # nunu final state
    #######
    df_nunu = df.Filter(select_nunu)
    hists.append(df_nunu.Histo1D(("cutFlow", "", *bins_count), "cut11"))
    hists.append(df_nunu.Histo1D(("znunu_haa_m", "", *bins_m_zoom), "haa_m"))

    # extra variables
    hists.append(df_nunu.Histo1D(("znunu_cosThetaMiss", "", *bins_cosThetaMiss), "cosThetaMiss"))
    hists.append(df_nunu.Histo1D(("znunu_missingEnergy", "", *bins_m), "missingEnergy"))
    hists.append(df_nunu.Histo1D(("znunu_cosThetaPhotons", "", *(100, 0, 1)), "cosThetaPhotons"))



    #######
    # mumu final state
    #######
    df_mumu = df.Filter(select_mumu)
    hists.append(df_mumu.Histo1D(("cutFlow", "", *bins_count), "cut12"))

    df_mumu = df_mumu.Define("muons_tlv", "FCCAnalyses::makeLorentzVectors(muons)")
    df_mumu = df_mumu.Define("dimuon", "muons_tlv[0]+muons_tlv[1]")
    df_mumu = df_mumu.Define("dimuon_m", "dimuon.M()")
    df_mumu = df_mumu.Define("dimuon_p", "dimuon.P()")
    df_mumu = df_mumu.Define("dimuon_theta", "dimuon.Theta()")


    hists.append(df_mumu.Histo1D(("zmumu_m_nOne", "", *bins_m), "dimuon_m"))
    df_mumu = df_mumu.Filter("dimuon_m > 80 && dimuon_m < 100")
    hists.append(df_mumu.Histo1D(("zmumu_p_nOne", "", *bins_m), "dimuon_p"))
    hists.append(df_mumu.Histo1D(("zmumu_haa_m", "", *bins_m_zoom), "haa_m"))


    df_mumu = df_mumu.Define("zmumu_dtheta", "dimuon.Theta() - haa_tlv.Theta()")
    df_mumu = df_mumu.Define("zmumu_dphi", "ROOT::Math::VectorUtil::DeltaPhi(dimuon, haa_tlv)")
    hists.append(df_mumu.Histo1D(("zmumu_dtheta", "", *(200, -10, 10)), "zmumu_dtheta"))
    hists.append(df_mumu.Histo1D(("zmumu_dphi", "", *(200, -10, 10)), "zmumu_dphi"))



    # extra variables

    # cos(theta1,2) CP analysis inspired
    df_mumu = df_mumu.Define("costheta1", "abs(cos(dimuon_theta))")
    df_mumu = df_mumu.Define("muon_neg_tlv", "(muons_q[0] < 0) ? muons_tlv[0] : muons_tlv[1]")
    df_mumu = df_mumu.Define("vec3_muon_neg", "ROOT::Math::XYZVector(muon_neg_tlv.Px(),muon_neg_tlv.Py(),muon_neg_tlv.Pz());")
    df_mumu = df_mumu.Define("vec3_mumu", "ROOT::Math::XYZVector(dimuon.Px(),dimuon.Py(),dimuon.Pz());")
    df_mumu = df_mumu.Define("costheta2", "abs(cos((vec3_muon_neg - vec3_mumu).Theta()))")
    hists.append(df_mumu.Histo1D(("zmumu_costheta1", "", *(100, -1, 1)), "costheta1"))
    hists.append(df_mumu.Histo1D(("zmumu_costheta2", "", *(100, -1, 1)), "costheta2"))

    # angle between Higgs and Z
    df_mumu = df_mumu.Define("photon_mumu_costheta", "abs(cos((vec3_diphoton-vec3_mumu).Theta()))")
    hists.append(df_mumu.Histo1D(("zmumu_photon_mumu_costheta", "", *(100, -1, 1)), "photon_mumu_costheta"))




    #######
    # ee final state
    #######
    df_ee = df.Filter(select_ee)
    hists.append(df_ee.Histo1D(("cutFlow", "", *bins_count), "cut13"))

    df_ee = df_ee.Define("electrons_tlv", "FCCAnalyses::makeLorentzVectors(electrons)")
    df_ee = df_ee.Define("dielectron", "electrons_tlv[0]+electrons_tlv[1]")
    df_ee = df_ee.Define("dielectron_m", "dielectron.M()")

    hists.append(df_ee.Histo1D(("zee_m_nOne", "", *bins_m), "dielectron_m"))
    df_ee = df_ee.Filter("dielectron_m > 80 && dielectron_m < 100")
    hists.append(df_ee.Histo1D(("zee_m", "", *bins_m), "dielectron_m"))
    hists.append(df_ee.Histo1D(("zee_haa_m", "", *bins_m_zoom), "haa_m"))


    return hists, weightsum
