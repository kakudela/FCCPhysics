
import ROOT
ROOT.TH1.SetDefaultSumw2(ROOT.kTRUE)



# list of all processes
fraction = 1
processList = {
    'wzp6_ee_nunuH_HZa_ecm240':     {'fraction':1},
    'wzp6_ee_eeH_HZa_ecm240':       {'fraction':1},
    'wzp6_ee_tautauH_HZa_ecm240':   {'fraction':1},
    'wzp6_ee_ccH_HZa_ecm240':       {'fraction':1},
    'wzp6_ee_bbH_HZa_ecm240':       {'fraction':1},
    'wzp6_ee_qqH_HZa_ecm240':       {'fraction':1},
    'wzp6_ee_ssH_HZa_ecm240':       {'fraction':1},
    'wzp6_ee_mumuH_HZa_ecm240':     {'fraction':1},
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
    #'wz3p6_ee_ee_Mee_30p_ecm240':   {'fraction':fraction},
    'wz3p6_ee_nunu_ecm240':         {'fraction':fraction},
}





inputDir = "/ceph/submit/data/group/fcc/ee/generation/DelphesEvents/winter2023/IDEA/"
procDict = "/ceph/submit/data/group/fcc/ee/generation/DelphesEvents/winter2023/IDEA/samplesDict.json"

# additional/custom C++ functions
includePaths = ["../../functions/functions.h", "../../functions/functions_gen.h"]


# output directory
outputDir   = "output/h_za/histmaker/"


# optional: ncpus, default is 4, -1 uses all cores available
nCPUS       = 32

# scale the histograms with the cross-section and integrated luminosity
doScale = True
intLumi = 10800000

# define histograms
bins_m = (250, 0, 250) # 100 MeV bins
bins_p = (200, 0, 200) # 100 MeV bins
bins_m_zoom = (200, 110, 130) # 100 MeV

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

bins_cos = (20000, -1, 1)
bins_cos_abs = (10000, 0, 1)
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

    # get all photons
    df = df.Alias("Photon0", "Photon#0.index")
    df = df.Define("photons_all", "FCCAnalyses::ReconstructedParticle::get(Photon0, ReconstructedParticles)")
    df = df.Define("photons_all_p", "FCCAnalyses::ReconstructedParticle::get_p(photons_all)")
    df = df.Define("photons_all_n", "FCCAnalyses::ReconstructedParticle::get_n(photons_all)")

    df_ph = df.Filter("photons_all_n > 0")
    df_ph = df_ph.Define("leading_photon_p", "photons_all_p[0]")
    hists.append(df_ph.Histo1D(("leading_photon_p", "", *bins_m), "leading_photon_p"))



    df = df.Define("photons", "FCCAnalyses::sel_range(20, 50, false)(photons_all, photons_all_p)") # based on vvH sample
    df = df.Define("photons_p", "FCCAnalyses::ReconstructedParticle::get_p(photons)")
    df = df.Define("photons_theta", "FCCAnalyses::ReconstructedParticle::get_theta(photons)")
    df = df.Define("photons_n", "FCCAnalyses::ReconstructedParticle::get_n(photons)")
    df = df.Define("photons_costheta", "FCCAnalyses::Vec_f ret; for(auto & theta: photons_theta) ret.push_back(std::abs(cos(theta))); return ret;")
    df = df.Define("photons_iso", "FCCAnalyses::coneIsolation(0.01, 0.5)(photons, ReconstructedParticles)")


    df = df.Define("photon", "photons[0]")
    df = df.Define("photon_vec", "Vec_rp{photon}")
    df = df.Define("photon_tlv", "FCCAnalyses::makeLorentzVectors(photon_vec)[0]")
    df = df.Define("photon_p", "photons_p[0]")
    df = df.Define("photon_iso", "photons_iso[0]")
    df = df.Define("photon_costheta", "photons_costheta[0]")


    #########
    ### CUT 0: all events
    #########
    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut0"))


    #########
    ### CUT 1: at least 1 photon
    #########
    hists.append(df.Histo1D(("photons_n", "", *bins_count), "photons_n"))
    df = df.Filter("photons_n >= 1")
    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut1"))


    #########
    ### CUT 2: exactly 1 photon
    #########
    df = df.Filter("photons_n == 1")
    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut2"))


    #########
    ### CUT 3: photon isolation
    #########
    hists.append(df.Histo1D(("photon_iso", "", *bins_iso), "photon_iso"))
    df = df.Filter("photon_iso < 0.3")
    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut3"))


    #########
    ### CUT 4: photon cos(theta) 
    #########
    hists.append(df.Histo1D(("photon_costheta", "", *bins_cos_abs), "photon_costheta"))
    df = df.Filter("photon_costheta < 0.85")
    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut4"))


    ####
    ## CUT 5: cos theta(miss)
    ####
    df = df.Define("missingEnergy_rp", "FCCAnalyses::missingEnergy(240., ReconstructedParticles)")
    df = df.Define("missingEnergy_tlv", "FCCAnalyses::makeLorentzVectors(missingEnergy_rp)[0]")
    df = df.Define("missingEnergy", "missingEnergy_rp[0].energy")
    df = df.Define("cosThetaEmiss", "FCCAnalyses::get_cosTheta_miss(missingEnergy_rp)")
    hists.append(df.Histo1D(("cosThetaEmiss", "", *bins_cosThetaMiss), "cosThetaEmiss"))
    df = df.Filter("cosThetaEmiss < 0.995")
    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut5"))

    df = df.Define("missingMass", "FCCAnalyses::missingMass(240., ReconstructedParticles)")
    hists.append(df.Histo1D(("missingEnergy", "", *bins_m), "missingEnergy"))
    hists.append(df.Histo1D(("missingMass", "", *bins_m), "missingMass"))

    df = df.Define("rp_no", "ReconstructedParticles.size()")
    hists.append(df.Histo1D(("rp_no", "", *(200, 0, 200)), "rp_no"))
    



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
    hists.append(df.Histo1D(("muons_all_p_cut0", "", *bins_p), "muons_all_p"))
    hists.append(df.Histo1D(("muons_all_theta_cut0", "", *bins_theta), "muons_all_theta"))
    hists.append(df.Histo1D(("muons_all_phi_cut0", "", *bins_phi), "muons_all_phi"))
    hists.append(df.Histo1D(("muons_all_q_cut0", "", *bins_charge), "muons_all_q"))
    hists.append(df.Histo1D(("muons_all_no_cut0", "", *bins_count), "muons_all_no"))

    hists.append(df.Histo1D(("electrons_all_p_cut0", "", *bins_p), "electrons_all_p"))
    hists.append(df.Histo1D(("electrons_all_theta_cut0", "", *bins_theta), "electrons_all_theta"))
    hists.append(df.Histo1D(("electrons_all_phi_cut0", "", *bins_phi), "electrons_all_phi"))
    hists.append(df.Histo1D(("electrons_all_q_cut0", "", *bins_charge), "electrons_all_q"))
    hists.append(df.Histo1D(("electrons_all_no_cut0", "", *bins_count), "electrons_all_no"))

    #hists.append(df.Histo1D(("photons_p", "", *bins_p), "photons_p"))
    #hists.append(df.Histo1D(("photons_n", "", *bins_count), "photons_n"))


    
    


    ##### CATEGORIZATION: based on #muons, # electrons, missing energy
    #select_vv_mumu = "muons_no == 2 && electrons_no == 0 && (missingEnergy > 90 && missingEnergy < 180)" # cuts based on nunuHZa bare
    #select_vv_ee = "electrons_no == 2 && muons_no == 0 && (missingEnergy > 90 && missingEnergy < 180)"
    #select_vv_qq = "electrons_no == 0 && muons_no == 0 && (missingEnergy > 90 && missingEnergy < 180)"
    #select_vv_vv = "electrons_no == 0 && muons_no == 0 && (missingEnergy > 180)"
    
    
    select_vv_mumu = "muons_no == 2 && electrons_no == 0 && (missingEnergy > 90 && missingEnergy < 180)" # cuts based on nunuHZa bare
    select_vv_ee = "electrons_no == 2 && muons_no == 0 && (missingEnergy > 90 && missingEnergy < 180)"
    select_vv_qq = "electrons_no == 0 && muons_no == 0 && (missingEnergy > 90 && missingEnergy < 180)"
    select_vv_vv = "electrons_no == 0 && muons_no == 0 && (missingEnergy > 180)"

    select_vv_qq = "electrons_no == 0 && muons_no == 0 && (missingMass > (91.2-20) && missingMass < (91.2+20))"



    #######
    # vv_qq final state
    #######
    df_qq = df.Filter(select_vv_qq)
    hists.append(df_qq.Histo1D(("cutFlow", "", *bins_count), "cut6"))
    hists.append(df_qq.Histo1D(("vvqq_rp_no", "", *(200, 0, 200)), "rp_no"))
    
    df_qq = df_qq.Filter("rp_no > 15") # kills further nunu and tautau. Makes clustering healthy and well-defined
    
    #hists.append(df_qq.Histo1D(("zqq_photons_p", "", *bins_p), "photons_p"))
    #hists.append(df_qq.Histo1D(("zqq_photon_leading_costheta", "", *bins_cosThetaMiss), "photon_leading_costheta"))
    #hists.append(df_qq.Histo1D(("zqq_photons_costheta", "", *bins_cosThetaMiss), "photons_costheta"))
    #hists.append(df_qq.Histo1D(("zqq_photon_leading_p", "", *bins_p), "photon_leading_p"))

    # define PF candidates collection by removing the muons
    df_qq = df_qq.Define("rps_no_photon", "FCCAnalyses::ReconstructedParticle::remove(ReconstructedParticles, photon_vec)")
    df_qq = df_qq.Define("RP_px", "FCCAnalyses::ReconstructedParticle::get_px(rps_no_photon)")
    df_qq = df_qq.Define("RP_py", "FCCAnalyses::ReconstructedParticle::get_py(rps_no_photon)")
    df_qq = df_qq.Define("RP_pz","FCCAnalyses::ReconstructedParticle::get_pz(rps_no_photon)")
    df_qq = df_qq.Define("RP_e", "FCCAnalyses::ReconstructedParticle::get_e(rps_no_photon)")
    df_qq = df_qq.Define("RP_m", "FCCAnalyses::ReconstructedParticle::get_mass(rps_no_photon)")
    df_qq = df_qq.Define("RP_q", "FCCAnalyses::ReconstructedParticle::get_charge(rps_no_photon)")
    df_qq = df_qq.Define("pseudo_jets", "FCCAnalyses::JetClusteringUtils::set_pseudoJets(RP_px, RP_py, RP_pz, RP_e)")

    df_qq = df_qq.Define("clustered_jets", "JetClustering::clustering_ee_kt(2, 2, 0, 10)(pseudo_jets)")
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
    df_qq = df_qq.Define("jet1_nc", "jetconstituents[0].size()")
    df_qq = df_qq.Define("jet2_nc", "jetconstituents[0].size()")
    df_qq = df_qq.Define("qq", "jet1+jet2")
    df_qq = df_qq.Define("qq_tlv", "TLorentzVector ret; ret.SetPxPyPzE(qq.Px(), qq.Py(), qq.Pz(), qq.E()); return ret;")
    df_qq = df_qq.Define("qq_m", "qq.M()")
    df_qq = df_qq.Define("qq_p", "qq.P()")
    df_qq = df_qq.Define("qq_costheta", "std::cos(qq.Theta())")

    hists.append(df_qq.Histo1D(("vvqq_jet1_p", "", *bins_m), "jet1_p"))
    hists.append(df_qq.Histo1D(("vvqq_jet2_p", "", *bins_m), "jet2_p"))
    
    hists.append(df_qq.Histo1D(("vvqq_jet1_nc", "", *(100, 0, 100)), "jet1_nc"))
    hists.append(df_qq.Histo1D(("vvqq_jet2_nc", "", *(100, 0, 100)), "jet2_nc"))
    
    # jet energy -- seems redundant with the other cuts, but keep it
    df_qq = df_qq.Filter("jet1_p > 10 && jet2_p > 10")
    
    # first select good Z candidate based on mass
    hists.append(df_qq.Histo1D(("vvqq_qq_m", "", *bins_m), "qq_m"))
    df_qq = df_qq.Filter("qq_m < (91.2+15) && qq_m > (91.2-15)")
    
    
    hists.append(df_qq.Histo1D(("vvqq_qq_costheta", "", *(100, 0 , 1)), "qq_costheta"))
    df_qq = df_qq.Filter("qq_costheta < 0.7") # quite tight cut
    
    # then select whether the Z comes from the higgs or not


    # za system
    df_qq = df_qq.Define("za", "qq_tlv+photon_tlv")
    df_qq = df_qq.Define("za_m", "za.M()") # nice peak at 125 GeV for nunuH --> correct clustering and selection of leading photon
    df_qq = df_qq.Define("za_p", "za.P()")

    # recoil of za system
    df_qq = df_qq.Define("init_tlv", "ROOT::Math::PxPyPzEVector(0, 0, 0, 240)")
    df_qq = df_qq.Define("init_tlvv", "TLorentzVector ret; ret.SetPxPyPzE(init_tlv.Px(), init_tlv.Py(), init_tlv.Pz(), init_tlv.E()); return ret;")
    df_qq = df_qq.Define("recoil_za", "init_tlvv-za")
    df_qq = df_qq.Define("recoil_za_m", "recoil_za.M()")
    df_qq = df_qq.Define("recoil_za_p", "recoil_za.P()")

    df_qq = df_qq.Define("recoil_z", "init_tlvv-qq_tlv") # higgs
    df_qq = df_qq.Define("recoil_z_m", "recoil_z.M()")
    
    
    # vva system
    
    df_qq = df_qq.Define("missingMass_vec", "ROOT::Math::PxPyPzMVector(missingEnergy_tlv.Px(), missingEnergy_tlv.Py(), missingEnergy_tlv.Pz(), missingMass)")
    df_qq = df_qq.Define("missingMass_tlv", "TLorentzVector ret; ret.SetPxPyPzE(missingMass_vec.Px(), missingMass_vec.Py(), missingMass_vec.Pz(), missingMass_vec.E()); return ret;")
    df_qq = df_qq.Define("vva", "missingMass_tlv+photon_tlv")
    df_qq = df_qq.Define("vva_m", "vva.M()")
    df_qq = df_qq.Define("vva_p", "vva.P()")
    #df_qq = df_qq.Define("vva_m", "missingMass + photon_tlv.E()")

    # recoil of vva system
    df_qq = df_qq.Define("recoil_vva", "init_tlvv-vva")
    df_qq = df_qq.Define("recoil_vva_m", "recoil_vva.M()")
    df_qq = df_qq.Define("recoil_vva_p", "recoil_vva.P()")

    df_qq = df_qq.Define("recoil_vv", "init_tlvv-missingMass_tlv") # higgs
    df_qq = df_qq.Define("recoil_vv_m", "recoil_vv.M()")

    
    
    ## mass differences
    df_qq = df_qq.Define("mass_difference_qq", "za_m - qq_m")
    df_qq = df_qq.Define("mass_difference_vv", "vva_m - qq_m")
    
    
    # recoil_za_m peaks the same for vvqqa or qqvva --> generic cut
    hists.append(df_qq.Histo1D(("vvqq_recoil_za_m", "", *bins_m), "recoil_za_m"))
    hists.append(df_qq.Histo1D(("vvqq_recoil_vva_m", "", *bins_m), "recoil_vva_m"))
    df_qq = df_qq.Filter("recoil_za_m > 80 && recoil_za_m < 100")
    
    hists.append(df_qq.Histo1D(("vvqq_mass_difference_qq", "", *bins_m), "mass_difference_qq"))
    hists.append(df_qq.Histo1D(("vvqq_mass_difference_vv", "", *bins_m), "mass_difference_vv"))
    
    # now we need to split the events based on vvqqa or qqvva --> mass separation
    hists.append(df_qq.Histo1D(("vvqq_qq_p", "", *bins_m), "qq_p"))
    hists.append(df_qq.Histo1D(("vvqq_za_m", "", *bins_m), "za_m"))
    hists.append(df_qq.Histo1D(("vvqq_za_p", "", *bins_m), "za_p"))
    hists.append(df_qq.Histo1D(("vvqq_recoil_z_m", "", *bins_m), "recoil_z_m"))

    hists.append(df_qq.Histo1D(("vvqq_vva_m", "", *bins_m), "vva_m"))
    hists.append(df_qq.Histo1D(("vvqq_vva_p", "", *bins_m), "vva_p"))
    hists.append(df_qq.Histo1D(("vvqq_recoil_vv_m", "", *bins_m), "recoil_vv_m"))

    
    
    #df_qq = df_qq.Define("chi2", "std::pow(za_m-125, 2) + 0.0*std::pow(za_p-50, 2) - 0.0*std::pow(qq_p-50, 2) - std::pow(recoil_z_m-125, 2)")
    df_qq = df_qq.Define("chi2", "std::pow(za_m-125, 2) + 0.0*std::pow(recoil_vv_m-125, 2) - std::pow(recoil_z_m-125, 2) - std::pow(vva_m-125, 2)")
    # adding qq_p and za_p makes it worse
    hists.append(df_qq.Histo1D(("vvqq_chi2", "", *(2400, -200, 200)), "chi2"))
    
    hists.append(df_qq.Histo2D(("final_histo", "", *(30, 110, 140, 30, 110, 140)), "za_m", "vva_m"))
    
    # split based on vvH sample (== qqa)

    
    # split the events based on chi2 --> can also be done with classifier based on masses, momenta and angles?
    df_qq_qqa = df_qq.Filter("chi2 < 0") # qqa
    hists.append(df_qq_qqa.Histo1D(("vvqqqqa_za_m", "", *bins_m), "za_m")) # higgs
    hists.append(df_qq_qqa.Histo1D(("vvqqqqa_za_p", "", *bins_m), "za_p"))
    hists.append(df_qq_qqa.Histo1D(("vvqqqqa_recoil_za_m", "", *bins_m), "recoil_za_m")) # z mass
    
    hists.append(df_qq_qqa.Histo1D(("vvqqqqa_final", "", *(30, 110, 140)), "za_m")) # higgs mass


    df_qq_vva = df_qq.Filter("chi2 > 0") # vva
    hists.append(df_qq_vva.Histo1D(("vvqqvva_z_p", "", *bins_m), "qq_p"))
    hists.append(df_qq_vva.Histo1D(("vvqqvva_recoil_z_m", "", *bins_m), "recoil_z_m")) # higgs mass
    hists.append(df_qq_vva.Histo1D(("vvqqvva_final", "", *(30, 110, 140)), "recoil_z_m")) # higgs mass

    
    #hists.append(df_qq.Histo1D(("recoil_qq_m", "", *bins_m), "recoil_qq_m"))
    #hists.append(df_qq.Histo1D(("recoil_qqa_m", "", *bins_m), "recoil_qqa_m"))

    #hists.append(df_qq.Histo1D(("zqq_m", "", *bins_m), "dijet_m"))
    #hists.append(df_qq.Histo1D(("zqq_haa_m", "", *bins_m_zoom), "haa_m"))
    #hists.append(df_qq.Histo1D(("zqq_aa_p", "", *bins_p), "haa_p"))
    #hists.append(df_qq.Histo1D(("zqq_qq_p", "", *bins_p), "dijet_p"))

    #hists.append(df_qq.Histo1D(("zqq_acoplanarity", "", *bins_aco), "acoplanarity"))
    #hists.append(df_qq.Histo1D(("zqq_acolinearity", "", *bins_aco), "acolinearity"))
    #hists.append(df_qq.Histo1D(("zqq_cosThetaMiss", "", *bins_cosThetaMiss), "cosTheta_miss"))






    #######
    # vv_vv final state
    #######
    df_nunu = df.Filter(select_vv_vv)
    hists.append(df_nunu.Histo1D(("cutFlow", "", *bins_count), "cut11"))

    hists.append(df_nunu.Histo1D(("vv_vv_missingEnergy_nOne", "", *bins_m), "missingEnergy"))
    hists.append(df_nunu.Histo1D(("vv_vv_cosThetaMiss_nOne", "", *bins_cosThetaMiss), "cosThetaEmiss"))

    #hists.append(df_nunu.Histo1D(("znunu_haa_m", "", *bins_m_zoom), "haa_m"))


    #######
    # vv_mumu final state
    #######
    df_mumu = df.Filter(select_vv_mumu)
    hists.append(df_mumu.Histo1D(("cutFlow", "", *bins_count), "cut12"))

    df_mumu = df_mumu.Define("muons_tlv", "FCCAnalyses::makeLorentzVectors(muons)")
    df_mumu = df_mumu.Define("dimuon", "muons_tlv[0]+muons_tlv[1]")
    df_mumu = df_mumu.Define("dimuon_m", "dimuon.M()")

    hists.append(df_mumu.Histo1D(("vv_mumu_mumu_m_nOne", "", *bins_m), "dimuon_m"))
    #df_mumu = df_mumu.Filter("dimuon_m > 80 && dimuon_m < 100")
    #hists.append(df_mumu.Histo1D(("vv_mumu_mumu_m", "", *bins_m), "dimuon_m"))
    #hists.append(df_mumu.Histo1D(("zmumu_haa_m", "", *bins_m_zoom), "haa_m"))



    # za system
    df_mumu = df_mumu.Define("za", "dimuon+photon_tlv")
    df_mumu = df_mumu.Define("za_m", "za.M()")
    df_mumu = df_mumu.Define("za_p", "za.P()")

    # recoil system of 
    df_mumu = df_mumu.Define("init_tlv", "ROOT::Math::PxPyPzEVector(0, 0, 0, 240)")
    df_mumu = df_mumu.Define("init_tlvv", "TLorentzVector ret; ret.SetPxPyPzE(init_tlv.Px(), init_tlv.Py(), init_tlv.Pz(), init_tlv.E()); return ret;")
    df_mumu = df_mumu.Define("recoil_za", "init_tlvv-za")
    df_mumu = df_mumu.Define("recoil_za_m", "recoil_za.M()")

    df_mumu = df_mumu.Define("recoil_z", "init_tlvv-dimuon") # higgs
    df_mumu = df_mumu.Define("recoil_z_m", "recoil_z.M()")
    
    # chi2
    df_mumu = df_mumu.Define("chi2", "std::pow(za_m-125, 2) + std::pow(recoil_za_m-91.2, 2) + std::pow(recoil_z_m-125, 2)")
    hists.append(df_mumu.Histo1D(("vv_mumu_chi2", "", *(300, 0, 300)), "chi2"))
    
    hists.append(df_mumu.Histo1D(("vv_mumu_recoil_za_m", "", *bins_m), "recoil_za_m"))
    

    df_qq = df_qq.Filter("recoil_za_m > 80 && recoil_za_m < 100")

    hists.append(df_mumu.Histo1D(("vv_mumu_za_m", "", *bins_m), "za_m"))
    hists.append(df_mumu.Histo1D(("vv_mumu_recoil_z_m", "", *bins_m), "recoil_z_m"))

    #######
    # ee final state
    #######
    df_ee = df.Filter(select_vv_ee)
    hists.append(df_ee.Histo1D(("cutFlow", "", *bins_count), "cut13"))

    df_ee = df_ee.Define("electrons_tlv", "FCCAnalyses::makeLorentzVectors(electrons)")
    df_ee = df_ee.Define("dielectron", "electrons_tlv[0]+electrons_tlv[1]")
    df_ee = df_ee.Define("dielectron_m", "dielectron.M()")

    hists.append(df_ee.Histo1D(("vv_ee_ee_m_nOne", "", *bins_m), "dielectron_m"))
    #df_ee = df_ee.Filter("dielectron_m > 80 && dielectron_m < 100")
    #hists.append(df_ee.Histo1D(("zee_m", "", *bins_m), "dielectron_m"))
    #hists.append(df_ee.Histo1D(("zee_haa_m", "", *bins_m_zoom), "haa_m"))











    #df = df.Filter("cosTheta_miss < .99") # 0.98
    #hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut7"))

    return hists, weightsum

    ####
    ## CUT 6: acolinearity
    ####
    #hists.append(df.Histo1D(("acolinearity_nOne", "", *bins_aco), "acolinearity"))
    #df = df.Filter("acolinearity > 0.1") # 0.98
    #hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut6"))



    #########
    ### CUT 3 :at least one resonance
    ### Exactly 2 photons
    #########

    '''
    # build the H resonance based on the available muons. Returns the best muon pair compatible with the H mass and recoil at 91 GeV
    # technically, it returns a ReconstructedParticleData object with index 0 the di-lepton system (H), index and 2 the leptons of the pair
    df = df.Define("hbuilder_result", "FCCAnalyses::resonanceBuilder_mass_recoil(125, 91.2, 0.5, 240, false)(photons, MCRecoAssociations0, MCRecoAssociations1, ReconstructedParticles, Particle, Particle0, Particle1)")
    
    df = df.Filter("hbuilder_result.size() > 0")
    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut3"))
    
    df = df.Define("haa", "ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>{hbuilder_result[0]}") # the H
    df = df.Define("haa_photons", "ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>{hbuilder_result[1],hbuilder_result[2]}") # the photons
    df = df.Define("haa_m", "FCCAnalyses::ReconstructedParticle::get_mass(haa)[0]")
    df = df.Define("haa_p", "FCCAnalyses::ReconstructedParticle::get_p(haa)[0]")
    df = df.Define("haa_recoil", "FCCAnalyses::ReconstructedParticle::recoilBuilder(240)(haa)")
    df = df.Define("haa_recoil_m", "FCCAnalyses::ReconstructedParticle::get_mass(haa_recoil)[0]")
    '''
    
    df = df.Filter("photons_n == 2")
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
    df = df.Filter("haa_recoil_m > 80 && haa_recoil_m < 120")
    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut4"))


    #####
    ### CUT 5: momentum
    #####
    hists.append(df.Histo1D(("haa_p_nOne", "", *bins_p), "haa_p"))
    df = df.Filter("haa_p > 20 && haa_p < 65")
    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut5"))


    df = df.Define("acoplanarity", "FCCAnalyses::acoplanarity(haa_photons)")
    df = df.Define("acolinearity", "FCCAnalyses::acolinearity(haa_photons)")
    hists.append(df.Histo1D(("acoplanarity_nOne", "", *bins_aco), "acoplanarity"))
    


    ####
    ## CUT 6: acolinearity
    ####
    hists.append(df.Histo1D(("acolinearity_nOne", "", *bins_aco), "acolinearity"))
    df = df.Filter("acolinearity > 0.1") # 0.98
    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut6"))

    ####
    ## CUT 7: cos theta(miss)
    ####
    df = df.Define("missingEnergy_rp", "FCCAnalyses::missingEnergy(240., ReconstructedParticles)")
    df = df.Define("missingEnergy", "missingEnergy_rp[0].energy")
    df = df.Define("cosTheta_miss", "FCCAnalyses::get_cosTheta_miss(missingEnergy_rp)")
    hists.append(df.Histo1D(("cosThetaMiss_nOne", "", *bins_cosThetaMiss), "cosTheta_miss"))
    df = df.Filter("cosTheta_miss < .99") # 0.98
    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut7"))


    ####
    ## CUT 8: missing energy
    ####
    hists.append(df.Histo1D(("missingEnergy_nOne", "", *bins_m), "missingEnergy"))
    df = df.Filter("missingEnergy < 115 && (missingEnergy > 95 || missingEnergy < 30 )")
    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut8"))


    #########
    ### CUT 9 :cut on Higgs mass
    #########
    hists.append(df.Histo1D(("haa_m_nOne", "", *bins_m), "haa_m"))
    df = df.Filter("haa_m > 110 && haa_m < 130")
    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut9"))



    #df = df.Define("acoplanarity", "FCCAnalyses::acoplanarity(haa_photons)")
    #df = df.Define("acolinearity", "FCCAnalyses::acolinearity(haa_photons)")
    hists.append(df.Histo1D(("acoplanarity", "", *bins_aco), "acoplanarity"))
    hists.append(df.Histo1D(("acolinearity", "", *bins_aco), "acolinearity"))
    hists.append(df.Histo1D(("haa_m", "", *bins_m_zoom), "haa_m"))


    ##### CATEGORIZATION: based on #muons, # electrons, missing energy
    select_mumu = "muons_no == 2 && electrons_no == 0 && missingEnergy < 30"
    select_ee = "electrons_no == 2 && muons_no == 0 && missingEnergy < 30"
    select_nunu = "electrons_no == 0 && muons_no == 0 && missingEnergy > 95"
    select_qq = "electrons_no == 0 && muons_no == 0 && missingEnergy < 30"







    return hists, weightsum
