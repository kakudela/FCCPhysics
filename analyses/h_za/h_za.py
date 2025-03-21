
import ROOT
import array
from addons.TMVAHelper.TMVAHelper import TMVAHelperXGB
ROOT.TH1.SetDefaultSumw2(ROOT.kTRUE)

ecm = 240

# list of all processes
fraction = 1

treemaker = False


if ecm == 240:
    processList = {
        'wzp6_ee_nunuH_HZa_ecm240':         {'fraction':1},
        'wzp6_ee_eeH_HZa_ecm240':           {'fraction':1},
        'wzp6_ee_tautauH_HZa_ecm240':       {'fraction':1},
        'wzp6_ee_mumuH_HZa_ecm240':         {'fraction':1},
        'wzp6_ee_ccH_HZa_ecm240':           {'fraction':1},
        'wzp6_ee_bbH_HZa_ecm240':           {'fraction':1},
        'wzp6_ee_qqH_HZa_ecm240':           {'fraction':1},
        'wzp6_ee_ssH_HZa_ecm240':           {'fraction':1},

        'wz3p6_ee_gammagamma_ecm240':       {'fraction':fraction},
        'wz3p6_ee_uu_ecm240':               {'fraction':fraction},
        'wz3p6_ee_dd_ecm240':               {'fraction':fraction},
        'wz3p6_ee_cc_ecm240':               {'fraction':fraction},
        'wz3p6_ee_ss_ecm240':               {'fraction':fraction},
        'wz3p6_ee_bb_ecm240':               {'fraction':fraction},
        'wz3p6_ee_tautau_ecm240':           {'fraction':fraction},
        'wz3p6_ee_mumu_ecm240':             {'fraction':fraction},
        'wz3p6_ee_ee_Mee_30_150_ecm240':    {'fraction':fraction},
        'wz3p6_ee_nunu_ecm240':             {'fraction':fraction},

        'p8_ee_ZZ_ecm240':             {'fraction':1},
    }

    if treemaker:
        processList = {
            'wzp6_ee_nunuH_HZa_ecm240':         {'fraction':1},
            'wzp6_ee_ccH_HZa_ecm240':           {'fraction':1},
            'wzp6_ee_bbH_HZa_ecm240':           {'fraction':1},
            'wzp6_ee_qqH_HZa_ecm240':           {'fraction':1},
            'wzp6_ee_ssH_HZa_ecm240':           {'fraction':1},
            'p8_ee_ZZ_ecm240':                  {'fraction':1},
            #'wz3p6_ee_gammagamma_ecm240':       {'fraction':fraction},
            #'wz3p6_ee_uu_ecm240':               {'fraction':fraction},
            #'wz3p6_ee_dd_ecm240':               {'fraction':fraction},
            #'wz3p6_ee_cc_ecm240':               {'fraction':fraction},
            #'wz3p6_ee_ss_ecm240':               {'fraction':fraction},
            #'wz3p6_ee_bb_ecm240':               {'fraction':fraction},
        }

if ecm == 365:
    processList = {
        'wzp6_ee_nunuH_HZa_ecm365':         {'fraction':1},
        'wzp6_ee_eeH_HZa_ecm365':           {'fraction':1},
        'wzp6_ee_tautauH_HZa_ecm365':       {'fraction':1},
        'wzp6_ee_mumuH_HZa_ecm365':         {'fraction':1},
        'wzp6_ee_ccH_HZa_ecm365':           {'fraction':1},
        'wzp6_ee_bbH_HZa_ecm365':           {'fraction':1},
        'wzp6_ee_qqH_HZa_ecm365':           {'fraction':1},
        'wzp6_ee_ssH_HZa_ecm365':           {'fraction':1},

        'wz3p6_ee_gammagamma_ecm365':       {'fraction':fraction},
        'wz3p6_ee_uu_ecm365':               {'fraction':fraction},
        'wz3p6_ee_dd_ecm365':               {'fraction':fraction},
        'wz3p6_ee_cc_ecm365':               {'fraction':fraction},
        'wz3p6_ee_ss_ecm365':               {'fraction':fraction},
        'wz3p6_ee_bb_ecm365':               {'fraction':fraction},
        'wz3p6_ee_tautau_ecm365':           {'fraction':fraction},
        'wz3p6_ee_mumu_ecm365':             {'fraction':fraction},
        'wz3p6_ee_ee_Mee_30_150_ecm365':    {'fraction':fraction},
        'wz3p6_ee_nunu_ecm365':             {'fraction':fraction},
    }





inputDir = "/ceph/submit/data/group/fcc/ee/generation/DelphesEvents/winter2023/IDEA/"
procDict = "/ceph/submit/data/group/fcc/ee/generation/DelphesEvents/winter2023/IDEA/samplesDict.json"

# additional/custom C++ functions
includePaths = ["../../functions/functions.h", "../../functions/functions_gen.h", "utils.h"]


# output directory
if treemaker:
    outputDir   = f"/ceph/submit/data/group/fcc/ee/analyses/h_za/treemaker/ecm{ecm}/qqvv_pairing_chi2" 
else:
    outputDir   = f"output/h_za/histmaker/ecm{ecm}/"


# optional: ncpus, default is 4, -1 uses all cores available
nCPUS       = 64

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

bins_m_mass_difference = (150, -50, 150)

bins_cos = (20000, -1, 1)
bins_cos_abs = (10000, 0, 1)
bins_iso = (1000, 0, 10)
bins_aco = (1000,0,5)
bins_dr = (500,0,5)
bins_cosThetaMiss = (10000, 0, 1)

bins_m_fine = (500, 110, 130) # 100 MeV bins
bins_merge = (50000, 0, 50000)

bins_higgs_m = (50, 100, 150)

ROOT.EnableImplicitMT(nCPUS) # hack
tmva_helper = TMVAHelperXGB("/ceph/submit/data/group/fcc/ee/analyses/za/training/bdt_model_ecm240.root", "bdt_model")
tmva_helper_hqqa = TMVAHelperXGB("/ceph/submit/data/group/fcc/ee/analyses/za/training/bdt_model_hqqa_ecm240.root", "bdt_model_hqqa")
tmva_helper_hvva = TMVAHelperXGB("/ceph/submit/data/group/fcc/ee/analyses/za/training/bdt_model_hvva_ecm240.root", "bdt_model_hvva")
tmva_helper_qqvv = TMVAHelperXGB("/ceph/submit/data/group/fcc/ee/analyses/za/training/bdt_model_qqvv_ecm240.root", "bdt_model_qqvv")


tmva_helper_qqvv_pairing_chi2 = TMVAHelperXGB("/ceph/submit/data/group/fcc/ee/analyses/za/training/qqvv_pairing_chi2.root", "qqvv_pairing_chi2")
tmva_helper_qqvv_qqa_splitting_mva = TMVAHelperXGB("/ceph/submit/data/group/fcc/ee/analyses/za/training/qqvv_qqa_splitting_mva.root", "qqvv_qqa_splitting_mva")
tmva_helper_qqvv_vva_splitting_mva = TMVAHelperXGB("/ceph/submit/data/group/fcc/ee/analyses/za/training/qqvv_vva_splitting_mva.root", "qqvv_vva_splitting_mva")


def build_graph_za(df, dataset):

    hists, cols = [], []

    df = df.Define("weight", "1.0")
    df = df.Define("ecm", "240" if ecm == 240 else "365")
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
    df = df.Define("cut15", "15")
    df = df.Define("cut16", "16")


    # define collections
    df = df.Alias("Particle0", "Particle#0.index")
    df = df.Alias("Particle1", "Particle#1.index")
    df = df.Alias("MCRecoAssociations0", "MCRecoAssociations#0.index")
    df = df.Alias("MCRecoAssociations1", "MCRecoAssociations#1.index")

    df = df.Define("RP_px_all", "FCCAnalyses::ReconstructedParticle::get_px(ReconstructedParticles)")
    df = df.Define("RP_py_all", "FCCAnalyses::ReconstructedParticle::get_py(ReconstructedParticles)")
    df = df.Define("RP_pz_all","FCCAnalyses::ReconstructedParticle::get_pz(ReconstructedParticles)")

    # get all photons
    df = df.Alias("Photon0", "Photon#0.index")
    df = df.Define("photons_all", "FCCAnalyses::ReconstructedParticle::get(Photon0, ReconstructedParticles)")
    df = df.Define("photons_all_p", "FCCAnalyses::ReconstructedParticle::get_p(photons_all)")
    df = df.Define("photons_all_n", "FCCAnalyses::ReconstructedParticle::get_n(photons_all)")
    df = df.Define("photons_all_theta", "FCCAnalyses::ReconstructedParticle::get_theta(photons_all)")
    df = df.Define("photons_all_costheta", "FCCAnalyses::Vec_f ret; for(auto & theta: photons_all_theta) ret.push_back(std::abs(cos(theta))); return ret;")

    hists.append(df.Histo1D(("photons_all_p", "", *bins_m), "photons_all_p"))
    hists.append(df.Histo1D(("photons_all_n", "", *bins_count), "photons_all_n"))
    hists.append(df.Histo1D(("photons_all_costheta", "", *bins_cos_abs), "photons_all_costheta"))


    # plot leading photon to select the momentum range
    df_ph = df.Filter("photons_all_n > 0")
    df_ph = df_ph.Define("leading_photon_p", "photons_all_p[0]")
    hists.append(df_ph.Histo1D(("leading_photon_p_nOne", "", *bins_m), "leading_photon_p"))


    df = df.Define("photons", "FCCAnalyses::sel_range(20, 50, false)(photons_all, photons_all_p)") # based on vvH sample # 15-60 GeV gives worse results
    df = df.Define("photons_p", "FCCAnalyses::ReconstructedParticle::get_p(photons)")
    df = df.Define("photons_theta", "FCCAnalyses::ReconstructedParticle::get_theta(photons)")
    df = df.Define("photons_n", "FCCAnalyses::ReconstructedParticle::get_n(photons)")
    df = df.Define("photons_costheta", "FCCAnalyses::Vec_f ret; for(auto & theta: photons_theta) ret.push_back(std::abs(cos(theta))); return ret;")
    df = df.Define("photons_iso", "FCCAnalyses::coneIsolation(0.01, 0.5)(photons, ReconstructedParticles)")

    hists.append(df.Histo1D(("photons_p", "", *bins_m), "photons_p"))
    hists.append(df.Histo1D(("photons_n", "", *bins_count), "photons_n"))
    hists.append(df.Histo1D(("photons_costheta", "", *bins_cos_abs), "photons_costheta"))

    df = df.Define("photon", "photons[0]")
    df = df.Define("photon_vec", "Vec_rp{photon}")
    df = df.Define("photon_tlv", "FCCAnalyses::makeLorentzVectors(photon_vec)[0]")
    df = df.Define("photon_p", "photons_p[0]")
    df = df.Define("photon_theta", "photons_theta[0]")
    df = df.Define("photon_iso", "photons_iso[0]")
    df = df.Define("photon_costheta", "photons_costheta[0]")



    #########
    ### CUT 0: all events
    #########
    hists.append(df.Histo1D(("qqvv_cutFlow", "", *bins_count), "cut0"))
    hists.append(df.Histo1D(("qqqq_cutFlow", "", *bins_count), "cut0"))


    #########
    ### CUT 1: at least 1 hard photon
    #########
    hists.append(df.Histo1D(("photons_n_nOne", "", *bins_count), "photons_n"))
    df = df.Filter("photons_n >= 1")
    hists.append(df.Histo1D(("qqvv_cutFlow", "", *bins_count), "cut1"))
    hists.append(df.Histo1D(("qqqq_cutFlow", "", *bins_count), "cut1"))


    #########
    ### CUT 2: exactly 1 hard photon
    #########
    df = df.Filter("photons_n == 1")
    hists.append(df.Histo1D(("qqvv_cutFlow", "", *bins_count), "cut2"))
    hists.append(df.Histo1D(("qqqq_cutFlow", "", *bins_count), "cut2"))

    # photon veto
    df = df.Define("photons_veto_all", "FCCAnalyses::sel_range(5, 9999, false)(photons_all, photons_all_p)")
    df = df.Define("photons_veto", "FCCAnalyses::ReconstructedParticle::remove(photons_veto_all, photon_vec)") # remove the selected photon
    df = df.Define("photons_veto_iso", "FCCAnalyses::coneIsolation(0.01, 0.5)(photons_veto, ReconstructedParticles)")
    df = df.Define("photons_veto_iso_sel", "FCCAnalyses::sel_range(0, 0.15, false)(photons_veto, photons_veto_iso)")
    df = df.Define("photons_veto_p", "FCCAnalyses::ReconstructedParticle::get_p(photons_veto_iso_sel)")
    df = df.Define("photons_veto_theta", "FCCAnalyses::ReconstructedParticle::get_theta(photons_veto_iso_sel)")
    df = df.Define("photons_veto_n", "FCCAnalyses::ReconstructedParticle::get_n(photons_veto_iso_sel)")
    df = df.Define("photons_veto_costheta", "FCCAnalyses::Vec_f ret; for(auto & theta: photons_veto_theta) ret.push_back(std::abs(cos(theta))); return ret;")

    hists.append(df.Histo1D(("photons_veto_p", "", *bins_m), "photons_veto_p"))
    hists.append(df.Histo1D(("photons_veto_n", "", *bins_count), "photons_veto_n"))
    hists.append(df.Histo1D(("photons_veto_costheta", "", *bins_cos_abs), "photons_veto_costheta"))

    #df = df.Filter("photons_veto_n == 0") # doesn't help as events with a photon veto is < 1% (see photon veto multiplicity)

    #########
    ### CUT 3: photon isolation
    #########
    hists.append(df.Histo1D(("photon_iso_nOne", "", *bins_iso), "photon_iso"))
    #df = df.Filter("photon_iso < 0.15")
    df = df.Filter("photon_iso < 0.3")
    hists.append(df.Histo1D(("qqvv_cutFlow", "", *bins_count), "cut3"))
    hists.append(df.Histo1D(("qqqq_cutFlow", "", *bins_count), "cut3"))


    #########
    ### CUT 4: photon cos(theta) 
    #########
    hists.append(df.Histo1D(("photon_costheta_nOne", "", *bins_cos_abs), "photon_costheta"))
    df = df.Filter("photon_costheta < 0.85")
    hists.append(df.Histo1D(("qqvv_cutFlow", "", *bins_count), "cut4"))
    hists.append(df.Histo1D(("qqqq_cutFlow", "", *bins_count), "cut4"))

    ####
    ## CUT 5: cos theta(miss)
    ####
    df = df.Define("missingEnergy_rp", "FCCAnalyses::missingEnergy(ecm, ReconstructedParticles)")
    df = df.Define("missingEnergy_tlv", "FCCAnalyses::makeLorentzVectors(missingEnergy_rp)[0]")
    df = df.Define("missingEnergy", "missingEnergy_rp[0].energy")
    df = df.Define("cosThetaEmiss", "FCCAnalyses::get_cosTheta_miss(missingEnergy_rp)")
    hists.append(df.Histo1D(("cosThetaEmiss_nOne", "", *bins_cosThetaMiss), "cosThetaEmiss"))
    #df = df.Filter("cosThetaEmiss < 0.995") # was before
    df = df.Filter("cosThetaEmiss < 0.98") # based on significance - difference betwee vvH and qqH
    hists.append(df.Histo1D(("qqvv_cutFlow", "", *bins_count), "cut5"))
    hists.append(df.Histo1D(("qqqq_cutFlow", "", *bins_count), "cut5"))



    # muons
    df = df.Alias("Muon0", "Muon#0.index")
    df = df.Define("muons_all", "FCCAnalyses::ReconstructedParticle::get(Muon0, ReconstructedParticles)")
    df = df.Define("muons_all_p", "FCCAnalyses::ReconstructedParticle::get_p(muons_all)")
    df = df.Define("muons_all_theta", "FCCAnalyses::ReconstructedParticle::get_theta(muons_all)")
    df = df.Define("muons_all_phi", "FCCAnalyses::ReconstructedParticle::get_phi(muons_all)")
    df = df.Define("muons_all_q", "FCCAnalyses::ReconstructedParticle::get_charge(muons_all)")
    df = df.Define("muons_all_no", "FCCAnalyses::ReconstructedParticle::get_n(muons_all)")

    df = df.Define("muons", "FCCAnalyses::ReconstructedParticle::sel_p(10)(muons_all)") ## was 2
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

    df = df.Define("electrons", "FCCAnalyses::ReconstructedParticle::sel_p(10)(electrons_all)") ## was 2
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


    ####
    ## CUT 6: lepton veto
    ####
    df = df.Filter("electrons_no == 0 && muons_no == 0")
    hists.append(df.Histo1D(("qqvv_cutFlow", "", *bins_count), "cut6"))
    hists.append(df.Histo1D(("qqqq_cutFlow", "", *bins_count), "cut6"))


    ####
    ## CUT 7: missing mass
    ####
    df = df.Define("missingMass", "FCCAnalyses::missingMass(ecm, ReconstructedParticles)")
    df = df.Define("missingMass_vec", "ROOT::Math::PxPyPzMVector(missingEnergy_tlv.Px(), missingEnergy_tlv.Py(), missingEnergy_tlv.Pz(), missingMass)")
    df = df.Define("missingMass_tlv", "TLorentzVector ret; ret.SetPxPyPzE(missingMass_vec.Px(), missingMass_vec.Py(), missingMass_vec.Pz(), missingMass_vec.E()); return ret;")
    hists.append(df.Histo1D(("missingMass_nOne", "", *bins_m), "missingMass"))


    df = df.Define("rps_no_photon", "FCCAnalyses::ReconstructedParticle::remove(ReconstructedParticles, photon_vec)")
    df = df.Define("init_tlv", f"TLorentzVector ret; ret.SetPxPyPzE(0, 0, 0, {ecm}); return ret;")

    df_qqqq = df.Filter("missingMass < 30")
    hists.append(df_qqqq.Histo1D(("qqqq_cutFlow", "", *bins_count), "cut7"))

    df_qqvv = df.Filter("missingMass > (91.2-25) && missingMass < (91.2+25)") # this kills the Z+hard ISR
    hists.append(df_qqvv.Histo1D(("qqvv_cutFlow", "", *bins_count), "cut7"))






    ############################################################################
    ##### qqvv final state
    ############################################################################

    hists.append(df_qqvv.Histo1D(("qqvv_missingEnergy_nOne", "", *bins_m), "missingEnergy"))


    df_qqvv = df_qqvv.Define("RP_px", "FCCAnalyses::ReconstructedParticle::get_px(rps_no_photon)")
    df_qqvv = df_qqvv.Define("RP_py", "FCCAnalyses::ReconstructedParticle::get_py(rps_no_photon)")
    df_qqvv = df_qqvv.Define("RP_pz","FCCAnalyses::ReconstructedParticle::get_pz(rps_no_photon)")
    df_qqvv = df_qqvv.Define("RP_e", "FCCAnalyses::ReconstructedParticle::get_e(rps_no_photon)")
    df_qqvv = df_qqvv.Define("RP_m", "FCCAnalyses::ReconstructedParticle::get_mass(rps_no_photon)")
    df_qqvv = df_qqvv.Define("RP_q", "FCCAnalyses::ReconstructedParticle::get_charge(rps_no_photon)")
    df_qqvv = df_qqvv.Define("pseudo_jets", "FCCAnalyses::JetClusteringUtils::set_pseudoJets(RP_px, RP_py, RP_pz, RP_e)")

    df_qqvv = df_qqvv.Define("clustered_jets", "JetClustering::clustering_ee_kt(2, 2, 0, 10)(pseudo_jets)")
    df_qqvv = df_qqvv.Define("jets", "FCCAnalyses::JetClusteringUtils::get_pseudoJets(clustered_jets)")
    df_qqvv = df_qqvv.Define("jetconstituents", "FCCAnalyses::JetClusteringUtils::get_constituents(clustered_jets)")
    df_qqvv = df_qqvv.Define("jets_n", "jetconstituents.size()")
    df_qqvv = df_qqvv.Define("jets_e", "FCCAnalyses::JetClusteringUtils::get_e(jets)")
    df_qqvv = df_qqvv.Define("jets_px", "FCCAnalyses::JetClusteringUtils::get_px(jets)")
    df_qqvv = df_qqvv.Define("jets_py", "FCCAnalyses::JetClusteringUtils::get_py(jets)")
    df_qqvv = df_qqvv.Define("jets_pz", "FCCAnalyses::JetClusteringUtils::get_pz(jets)")
    df_qqvv = df_qqvv.Define("jets_m", "FCCAnalyses::JetClusteringUtils::get_m(jets)")


    df_qqvv = df_qqvv.Define("dmerge_01", "FCCAnalyses::JetClusteringUtils::get_exclusive_dmerge(clustered_jets, 0)")
    df_qqvv = df_qqvv.Define("dmerge_12", "FCCAnalyses::JetClusteringUtils::get_exclusive_dmerge(clustered_jets, 1)")
    df_qqvv = df_qqvv.Define("dmerge_23", "FCCAnalyses::JetClusteringUtils::get_exclusive_dmerge(clustered_jets, 2)")
    hists.append(df_qqvv.Histo1D(("qqvv_dmerge_01", "", *bins_merge), "dmerge_01"))
    hists.append(df_qqvv.Histo1D(("qqvv_dmerge_12", "", *bins_merge), "dmerge_12"))
    hists.append(df_qqvv.Histo1D(("qqvv_dmerge_23", "", *bins_merge), "dmerge_23"))

    df_qqvv = df_qqvv.Filter("dmerge_12 > 1500")
    hists.append(df_qqvv.Histo1D(("qqvv_cutFlow", "", *bins_count), "cut8"))

    #### 
    ## CUT8: JET CLUSTERING
    ####
    # require good jet clustering
    hists.append(df_qqvv.Histo1D(("qqvv_jets_n_nOne", "", *bins_count), "jets_n"))
    df_qqvv = df_qqvv.Filter("jets_n == 2")
    hists.append(df_qqvv.Histo1D(("qqvv_cutFlow", "", *bins_count), "cut9"))

    df_qqvv = df_qqvv.Define("jet1", "ROOT::Math::PxPyPzEVector(jets_px[0], jets_py[0], jets_pz[0], jets_e[0])")
    df_qqvv = df_qqvv.Define("jet2", "ROOT::Math::PxPyPzEVector(jets_px[1], jets_py[1], jets_pz[1], jets_e[1])")
    df_qqvv = df_qqvv.Define("jet1_p", "jet1.P()")
    df_qqvv = df_qqvv.Define("jet2_p", "jet2.P()")
    df_qqvv = df_qqvv.Define("jet1_theta", "jet1.Theta()")
    df_qqvv = df_qqvv.Define("jet2_theta", "jet2.Theta()")
    df_qqvv = df_qqvv.Define("jet1_nc", "jetconstituents[0].size()")
    df_qqvv = df_qqvv.Define("jet2_nc", "jetconstituents[1].size()")
    df_qqvv = df_qqvv.Define("qq", "jet1+jet2")
    df_qqvv = df_qqvv.Define("qq_tlv", "TLorentzVector ret; ret.SetPxPyPzE(qq.Px(), qq.Py(), qq.Pz(), qq.E()); return ret;")
    df_qqvv = df_qqvv.Define("qq_m", "(float)qq.M()")
    df_qqvv = df_qqvv.Define("qq_p", "(float)qq.P()")
    df_qqvv = df_qqvv.Define("qq_costheta", "std::cos(qq.Theta())")

    hists.append(df_qqvv.Histo1D(("qqvv_jet1_p_nOne", "", *bins_m), "jet1_p"))
    hists.append(df_qqvv.Histo1D(("qqvv_jet2_p_nOne", "", *bins_m), "jet2_p"))
    
    hists.append(df_qqvv.Histo1D(("qqvv_jet1_nc_nOne", "", *(100, 0, 100)), "jet1_nc"))
    hists.append(df_qqvv.Histo1D(("qqvv_jet2_nc_nOne", "", *(100, 0, 100)), "jet2_nc"))


    ####
    ## CUT9: JET ID
    ####
    # jet constituents
    df_qqvv = df_qqvv.Filter("jet1_nc > 8 && jet2_nc > 8 && jet1_p > 15 && jet2_p > 15")
    hists.append(df_qqvv.Histo1D(("qqvv_cutFlow", "", *bins_count), "cut10"))

    #### 
    ## CUT10: qq mass
    ####
    # first select good Z candidate based on mass
    hists.append(df_qqvv.Histo1D(("qqvv_qq_m_nOne", "", *bins_m), "qq_m"))
    df_qqvv = df_qqvv.Filter("qq_m < (91.2+15) && qq_m > (91.2-15)")
    #df_qqvv = df_qqvv.Filter("qq_m < 105 && qq_m > 80.0") # mass candidate
    hists.append(df_qqvv.Histo1D(("qqvv_cutFlow", "", *bins_count), "cut11"))


    hists.append(df_qqvv.Histo1D(("qqvv_qq_costheta_nOne", "", *(100, 0 , 1)), "qq_costheta"))




    ## make possible combinations between qq and vv paired with a

    # qqa system
    df_qqvv = df_qqvv.Define("qqa", "qq_tlv+photon_tlv")
    df_qqvv = df_qqvv.Define("qqa_m", "(float)qqa.M()") # peaks at mH for vvH (sharp peak), broader peak at mH for qqH
    df_qqvv = df_qqvv.Define("qqa_p", "(float)qqa.P()")

    # recoil of qqa system
    df_qqvv = df_qqvv.Define("qqa_recoil", "init_tlv-qqa")
    df_qqvv = df_qqvv.Define("qqa_recoil_m", "(float)qqa_recoil.M()") # peaks at mZ for vvH, also for qqH (symmetric)
    df_qqvv = df_qqvv.Define("qqa_recoil_p", "(float)qqa_recoil.P()")

    # recoil of the qq system
    df_qqvv = df_qqvv.Define("qq_recoil", "init_tlv-qq_tlv") # peaks at mH for qqH, what for vvH?
    df_qqvv = df_qqvv.Define("qq_recoil_m", "(float)qq_recoil.M()")
    
    
    # vva system
    df_qqvv = df_qqvv.Define("vva", "missingMass_tlv+photon_tlv")
    df_qqvv = df_qqvv.Define("vva_m", "(float)vva.M()")
    df_qqvv = df_qqvv.Define("vva_p", "(float)vva.P()")

    # recoil of vva system
    df_qqvv = df_qqvv.Define("vva_recoil", "init_tlv-vva")
    df_qqvv = df_qqvv.Define("vva_recoil_m", "vva_recoil.M()")
    df_qqvv = df_qqvv.Define("vva_recoil_p", "vva_recoil.P()")

    # recoil of vv system
    df_qqvv = df_qqvv.Define("vv_m", "(float)missingMass_tlv.M()")
    df_qqvv = df_qqvv.Define("vv_recoil", "init_tlv-missingMass_tlv") # higgs
    df_qqvv = df_qqvv.Define("vv_recoil_m", "(float)vv_recoil.M()")


    #### 
    ## CUT11: cut on recoil of the qqa system (qqa_recoil_m peaks the same for vvqqa or qqvva --> generic cut)
    ####
    hists.append(df_qqvv.Histo1D(("qqvv_qqa_recoil_m_nOne", "", *bins_m), "qqa_recoil_m"))
    #df_qqvv = df_qqvv.Filter("qqa_recoil_m > 80 && qqa_recoil_m < 100") ## significance
    #df_qqvv = df_qqvv.Filter("qqa_recoil_m > 60 && qqa_recoil_m < 115") ## significance
    #hists.append(df_qqvv.Histo1D(("qqvv_cutFlow", "", *bins_count), "cut12"))

    hists.append(df_qqvv.Histo1D(("qqvv_vva_recoil_m_nOne", "", *bins_m), "vva_recoil_m"))


    # now we need to split the events based on vvqqa or qqvva --> mass separation
    hists.append(df_qqvv.Histo1D(("qqvv_qq_p", "", *bins_m), "qq_p"))
    hists.append(df_qqvv.Histo1D(("qqvv_qqa_m", "", *bins_m), "qqa_m"))
    hists.append(df_qqvv.Histo1D(("qqvv_qqa_p", "", *bins_m), "qqa_p"))
    hists.append(df_qqvv.Histo1D(("qqvv_qq_recoil_m", "", *bins_m), "qq_recoil_m"))

    hists.append(df_qqvv.Histo1D(("qqvv_vva_m", "", *bins_m), "vva_m"))
    hists.append(df_qqvv.Histo1D(("qqvv_vva_p", "", *bins_m), "vva_p"))
    hists.append(df_qqvv.Histo1D(("qqvv_vv_recoil_m", "", *bins_m), "vv_recoil_m"))

    # angular variables
    
    # angle between Z and H
    df_qqvv = df_qqvv.Define("dr_vv_qqa", "missingMass_tlv.DeltaR(qqa)")
    df_qqvv = df_qqvv.Define("dr_qq_vva", "qq_tlv.DeltaR(vva)")
    
    # angle between photon and vv, qq
    df_qqvv = df_qqvv.Define("dr_a_qq", "photon_tlv.DeltaR(qq_tlv)")
    df_qqvv = df_qqvv.Define("dr_a_vv", "photon_tlv.DeltaR(missingMass_tlv)")

    # mass differences
    df_qqvv = df_qqvv.Define("mass_difference_qqa_vv", "qqa_m - vv_m")
    df_qqvv = df_qqvv.Define("mass_difference_vva_qq", "vva_m - qq_m")
    hists.append(df_qqvv.Histo1D(("qqvv_mass_difference_qqa_vv", "", *bins_m_mass_difference), "mass_difference_qqa_vv"))
    hists.append(df_qqvv.Histo1D(("qqvv_mass_difference_vva_qq", "", *bins_m_mass_difference), "mass_difference_vva_qq"))

    # cosine angles
    df_qqvv = df_qqvv.Define("cos_qqa", "abs(cos(qqa.Theta()))")
    df_qqvv = df_qqvv.Define("cos_qq", "abs(cos(qq_tlv.Theta()))")
    df_qqvv = df_qqvv.Define("cos_vva", "abs(cos(vva.Theta()))")

    # transverse/longitudinal components of missing energy/qq
    df_qqvv = df_qqvv.Define("vv_trans", "missingMass_tlv.Pt()")
    df_qqvv = df_qqvv.Define("qq_trans", "qq_tlv.Pt()")
    df_qqvv = df_qqvv.Define("vv_long", "missingMass_tlv.Pz()")
    df_qqvv = df_qqvv.Define("qq_long", "qq_tlv.Pz()")
    
    hists.append(df_qqvv.Histo1D(("qqvv_vv_trans", "", *bins_m), "vv_trans"))
    hists.append(df_qqvv.Histo1D(("qqvv_qq_trans", "", *bins_m), "qq_trans"))
    hists.append(df_qqvv.Histo1D(("qqvv_vv_long", "", *bins_m), "vv_long"))
    hists.append(df_qqvv.Histo1D(("qqvv_qq_long", "", *bins_m), "qq_long"))



    ########################################
    ## try pairing of jets based on chi2
    ## not baseline analysis
    ########################################

    # pair the photon to the Z boson s
    df_qqvv = df_qqvv.Alias("Z1", "qq_tlv")
    df_qqvv = df_qqvv.Alias("Z2", "missingMass_tlv")
    
    df_qqvv = df_qqvv.Define("h_z_pairs", "FCCAnalyses::pairing_H_Z_a(ecm, Z1, Z2, photon_tlv)")
    df_qqvv = df_qqvv.Define("H", "h_z_pairs[0]") # Higgs
    df_qqvv = df_qqvv.Define("Z", "h_z_pairs[1]") # associated Z
    df_qqvv = df_qqvv.Define("ZH", "h_z_pairs[2]") # Z from Higgs

    df_qqvv = df_qqvv.Define("H_m", "(float)H.M()")
    df_qqvv = df_qqvv.Define("H_p", "(float)H.P()")
    df_qqvv = df_qqvv.Define("Z_m", "(float)Z.M()")
    df_qqvv = df_qqvv.Define("Z_p", "(float)Z.P()")
    df_qqvv = df_qqvv.Define("ZH_m", "(float)ZH.M()")
    df_qqvv = df_qqvv.Define("ZH_p", "(float)ZH.P()")



    # angle between Z and H
    df_qqvv = df_qqvv.Define("dr_H_Z", "(float)H.DeltaR(Z)")
    
    # angle between photon and Higgs/Z
    df_qqvv = df_qqvv.Define("dr_a_H", "(float)photon_tlv.DeltaR(H)")
    df_qqvv = df_qqvv.Define("dr_a_Z", "(float)photon_tlv.DeltaR(Z)")

    # mass differences
    df_qqvv = df_qqvv.Define("mass_difference_H_Z", "(float)(H_m - Z_m)")

    # cosine angles
    df_qqvv = df_qqvv.Define("cos_H", "(float)abs(cos(H.Theta()))")
    df_qqvv = df_qqvv.Define("cos_Z", "(float)abs(cos(Z.Theta()))")
    
    # recoil of the Higgs
    df_qqvv = df_qqvv.Define("H_recoil", "init_tlv-H")
    df_qqvv = df_qqvv.Define("H_recoil_m", "(float)H_recoil.M()") # peaks at mZ for vvH, also for qqH (symmetric)
    df_qqvv = df_qqvv.Define("H_recoil_p", "(float)H_recoil.P()")
    
    # recoil of the Higgs
    df_qqvv = df_qqvv.Define("Z_recoil", "init_tlv-Z")
    df_qqvv = df_qqvv.Define("Z_recoil_m", "(float)Z_recoil.M()") # peaks at mZ for vvH, also for qqH (symmetric)
    df_qqvv = df_qqvv.Define("Z_recoil_p", "(float)Z_recoil.P()")

    hists.append(df_qqvv.Histo1D(("qqvv_H_m", "", *bins_m), "H_m"))
    hists.append(df_qqvv.Histo1D(("qqvv_H_p", "", *bins_m), "H_p"))
    hists.append(df_qqvv.Histo1D(("qqvv_Z_m", "", *bins_m), "Z_m"))
    hists.append(df_qqvv.Histo1D(("qqvv_Z_p", "", *bins_m), "Z_p"))
    hists.append(df_qqvv.Histo1D(("qqvv_ZH_m", "", *bins_m), "ZH_m"))
    hists.append(df_qqvv.Histo1D(("qqvv_ZH_p", "", *bins_m), "ZH_p"))
    hists.append(df_qqvv.Histo1D(("qqvv_H_recoil_m", "", *bins_m), "H_recoil_m"))
    hists.append(df_qqvv.Histo1D(("qqvv_H_recoil_p", "", *bins_m), "H_recoil_p"))
    hists.append(df_qqvv.Histo1D(("qqvv_Z_recoil_m", "", *bins_m), "Z_recoil_m"))
    hists.append(df_qqvv.Histo1D(("qqvv_Z_recoil_p", "", *bins_m), "Z_recoil_p"))
    hists.append(df_qqvv.Histo1D(("qqvv_mass_difference_H_Z", "", *bins_m_mass_difference), "mass_difference_H_Z"))


    df_qqvv = tmva_helper_qqvv_pairing_chi2.run_inference(df_qqvv, col_name="qqvv_mva_score_split_chi2")
    hists.append(df_qqvv.Histo1D(("qqvv_mva_score_split_chi2", "", *(1000, 0, 1)), "qqvv_mva_score_split_chi2"))

    #varbins_mva_ = [0, 0.92, 1] # 0.75 = max significance
    #varbins_m_ = list(range(110, 141, 1))
    #varbins_mva = array.array('d', varbins_mva_)
    #varbins_m = array.array('d', varbins_m_)
    #model_qqa = ROOT.RDF.TH2DModel("qqvv_hqqa_final_2D_mva", "", len(varbins_m_)-1, varbins_m, len(varbins_mva_)-1, varbins_mva)
    #hists.append(df_qqvv_qqa_mva.Histo2D(model_qqa, "qqa_m", "hqqa_mva_score"))
    #df_qqvv = df_qqvv.Filter("qqvv_mva_score_split_chi2[0] > 0.75")
    hists.append(df_qqvv.Histo2D(("qqvv_qqa_m_vva_m", "", *((50, 100, 150)+(50, 100, 150))), "qqa_m", "vva_m"))



    ## variables for MVA
    mva_vars = ["H_m", "H_p", "Z_m", "Z_p", "ZH_m", "ZH_p", "jet1_theta", "jet2_theta", "jet1_p", "jet2_p", "photon_p", "photon_theta", "H_recoil_m", "H_recoil_p", "Z_recoil_m", "Z_recoil_p", "mass_difference_H_Z", "dr_H_Z", "dr_a_H", "dr_a_Z", "cos_H", "cos_Z", "vv_trans", "qq_trans", "vv_long", "qq_long"]
    return hists, weightsum, df_qqvv, mva_vars

    #return hists, weightsum, df, []

    ########################################
    ## split events based on chi2 or MVA whether the qq comes from the Higgs Z (hqq) or from the associated Z (zqq)
    ## chi2 not baseline analysis
    ########################################


    #df_qqvv = df_qqvv.Define("chi2", "std::pow(za_m-125, 2) + 0.0*std::pow(za_p-50, 2) - 0.0*std::pow(qq_p-50, 2) - std::pow(recoil_z_m-125, 2)")
    df_qqvv = df_qqvv.Define("chi2", "std::pow(qqa_m-125, 2) + std::pow(vv_recoil_m-125, 2) - std::pow(qq_recoil_m-125, 2) - std::pow(vva_m-125, 2)")
    # adding qq_p and za_p makes it worse
    hists.append(df_qqvv.Histo1D(("qqvv_chi2", "", *(2400, -200, 200)), "chi2"))



    # split the events based on chi2
    df_qqvv_qqa_chi2 = df_qqvv.Filter("chi2 < 0") # H(qqa)
    hists.append(df_qqvv_qqa_chi2.Histo1D(("qqvv_cutFlow", "", *bins_count), "cut14"))
    hists.append(df_qqvv_qqa_chi2.Histo1D(("qqvv_hqqa_qqa_m_chi2", "", *bins_m), "qqa_m")) # higgs mass
    hists.append(df_qqvv_qqa_chi2.Histo1D(("qqvv_hqqa_qqa_p_chi2", "", *bins_m), "qqa_p"))
    hists.append(df_qqvv_qqa_chi2.Histo1D(("qqvv_hqqa_qqa_recoil_m_chi2", "", *bins_m), "qqa_recoil_m")) # z mass
    hists.append(df_qqvv_qqa_chi2.Histo1D(("qqvv_hqqa_final_chi2", "", *(30, 110, 140)), "qqa_m")) # higgs mass
    hists.append(df_qqvv_qqa_chi2.Histo1D(("qqvv_hqqa_mass_difference_qqa_vv_chi2", "", *bins_m_mass_difference), "mass_difference_qqa_vv"))
    hists.append(df_qqvv_qqa_chi2.Histo1D(("qqvv_hqqa_mass_difference_vva_qq_chi2", "", *bins_m_mass_difference), "mass_difference_vva_qq"))

    df_qqvv_vva_chi2 = df_qqvv.Filter("chi2 > 0") # H(vva)
    hists.append(df_qqvv_vva_chi2.Histo1D(("qqvv_cutFlow", "", *bins_count), "cut15"))
    hists.append(df_qqvv_vva_chi2.Histo1D(("qqvv_hvva_qq_p_chi2", "", *bins_m), "qq_p"))
    hists.append(df_qqvv_vva_chi2.Histo1D(("qqvv_hvva_qq_recoil_m_chi2", "", *bins_m), "qq_recoil_m")) # higgs mass
    hists.append(df_qqvv_vva_chi2.Histo1D(("qqvv_hvva_final_chi2", "", *(30, 110, 140)), "qq_recoil_m")) # higgs mass
    hists.append(df_qqvv_vva_chi2.Histo1D(("qqvv_hvva_mass_difference_qqa_vv_chi2", "", *bins_m_mass_difference), "mass_difference_qqa_vv"))
    hists.append(df_qqvv_vva_chi2.Histo1D(("qqvv_hvva_mass_difference_vva_qq_chi2", "", *bins_m_mass_difference), "mass_difference_vva_qq"))


    ########################################
    ## split events based on MVA whether the qq comes from the Higgs Z (hqq) or from the associated Z (zqq)
    ## MVA not baseline analysis
    ########################################

    # split on MVA: signal=wzp6_ee_nunuH_HZa_ecm240, backgrounds = wzp6_ee_qqH_HZa_ecm240
    df_qqvv = tmva_helper.run_inference(df_qqvv, col_name="mva_score_qqa_vva")
    hists.append(df_qqvv.Histo1D(("qqvv_mva_score_qqa_vva", "", *(1000, 0, 1)), "mva_score_qqa_vva"))

    df_qqvv_qqa_mva = df_qqvv.Filter("mva_score_qqa_vva[0] > 0.5") # qqa
    hists.append(df_qqvv_qqa_mva.Histo1D(("qqvv_hqqa_final_1D_mva", "", *(30, 110, 140)), "qqa_m")) # higgs mass
    hists.append(df_qqvv_qqa_mva.Histo1D(("qqvv_cutFlow", "", *bins_count), "cut12"))

    hists.append(df_qqvv_qqa_mva.Histo1D(("qqvv_hqqa_qq_m_mva", "", *bins_m), "qq_m")) # z mass
    hists.append(df_qqvv_qqa_mva.Histo1D(("qqvv_hqqa_qq_p_mva", "", *bins_m), "qq_p"))
    hists.append(df_qqvv_qqa_mva.Histo1D(("qqvv_hqqa_qqa_m_mva", "", *bins_m), "qqa_m")) # higgs mass
    hists.append(df_qqvv_qqa_mva.Histo1D(("qqvv_hqqa_qqa_p_mva", "", *bins_m), "qqa_p"))
    hists.append(df_qqvv_qqa_mva.Histo1D(("qqvv_hqqa_qq_recoil_m_mva", "", *bins_m), "qq_recoil_m")) # ??
    hists.append(df_qqvv_qqa_mva.Histo1D(("qqvv_hqqa_qqa_recoil_za_m_mva", "", *bins_m), "qqa_recoil_m")) # z mass

    hists.append(df_qqvv_qqa_mva.Histo1D(("qqvv_hqqa_vv_m_mva", "", *bins_m), "vv_m"))
    hists.append(df_qqvv_qqa_mva.Histo1D(("qqvv_hqqa_vva_m_mva", "", *bins_m), "vva_m"))
    hists.append(df_qqvv_qqa_mva.Histo1D(("qqvv_hqqa_vva_p_mva", "", *bins_m), "vva_p"))
    hists.append(df_qqvv_qqa_mva.Histo1D(("qqvv_hqqa_vv_recoil_m_mva", "", *bins_m), "vv_recoil_m"))

    ## mass differences
    hists.append(df_qqvv_qqa_mva.Histo1D(("qqvv_hqqa_mass_difference_qqa_vv_mva", "", *bins_m_mass_difference), "mass_difference_qqa_vv"))
    hists.append(df_qqvv_qqa_mva.Histo1D(("qqvv_hqqa_mass_difference_vva_qq_mva", "", *bins_m_mass_difference), "mass_difference_vva_qq"))


    hists.append(df_qqvv_qqa_mva.Histo1D(("qqvv_hqqa_dr_vv_qqa_nOne", "", *bins_dr), "dr_vv_qqa"))
    hists.append(df_qqvv_qqa_mva.Histo1D(("qqvv_hqqa_dr_qq_vva_nOne", "", *bins_dr), "dr_qq_vva"))
    hists.append(df_qqvv_qqa_mva.Histo1D(("qqvv_hqqa_dr_a_qq_nOne", "", *bins_dr), "dr_a_qq"))
    hists.append(df_qqvv_qqa_mva.Histo1D(("qqvv_hqqa_dr_a_vv_nOne", "", *bins_dr), "dr_a_vv"))

    mva_vars = ["H_m", "H_p", "Z_m", "Z_p", "ZH_m", "ZH_p", "jet1_theta", "jet2_theta", "jet1_p", "jet2_p", "photon_p", "photon_theta", "H_recoil_m", "H_recoil_p", "Z_recoil_m", "Z_recoil_p", "mass_difference_H_Z", "dr_H_Z", "dr_a_H", "dr_a_Z", "cos_H", "cos_Z", "vv_trans", "qq_trans", "vv_long", "qq_long"]
    #return hists, weightsum, df_qqvv_qqa_mva, mva_vars


    df_qqvv_qqa_mva = tmva_helper_qqvv_qqa_splitting_mva.run_inference(df_qqvv_qqa_mva, col_name="hqqa_mva_score")
    hists.append(df_qqvv_qqa_mva.Histo1D(("qqvv_hqqa_mva_score", "", *(1000, 0, 1)), "hqqa_mva_score"))


    varbins_mva_ = [0, 0.92, 1] # 0.75 = max significance
    varbins_m_ = list(range(110, 141, 1))
    varbins_mva = array.array('d', varbins_mva_)
    varbins_m = array.array('d', varbins_m_)
    model_qqa = ROOT.RDF.TH2DModel("qqvv_hqqa_final_2D_mva", "", len(varbins_m_)-1, varbins_m, len(varbins_mva_)-1, varbins_mva)
    hists.append(df_qqvv_qqa_mva.Histo2D(model_qqa, "qqa_m", "hqqa_mva_score"))



    df_qqvv_vva_mva = df_qqvv.Filter("mva_score_qqa_vva[0] < 0.5") # vva
    hists.append(df_qqvv_vva_mva.Histo1D(("qqvv_hvva_final_1D_mva", "", *(30, 110, 140)), "qq_recoil_m")) # higgs mass
    hists.append(df_qqvv_vva_mva.Histo1D(("qqvv_cutFlow", "", *bins_count), "cut13"))

    hists.append(df_qqvv_vva_mva.Histo1D(("qqvv_hvva_qq_m_mva", "", *bins_m), "qq_m")) # z mass
    hists.append(df_qqvv_vva_mva.Histo1D(("qqvv_hvva_qq_p_mva", "", *bins_m), "qq_p"))
    hists.append(df_qqvv_vva_mva.Histo1D(("qqvv_hvva_qqa_m_mva", "", *bins_m), "qqa_m")) # higgs mass
    hists.append(df_qqvv_vva_mva.Histo1D(("qqvv_hvva_qqa_p_mva", "", *bins_m), "qqa_p"))
    hists.append(df_qqvv_vva_mva.Histo1D(("qqvv_hvva_qq_recoil_m_mva", "", *bins_m), "qq_recoil_m")) # ??
    hists.append(df_qqvv_vva_mva.Histo1D(("qqvv_hvva_qqa_recoil_za_m_mva", "", *bins_m), "qqa_recoil_m")) # z mass

    hists.append(df_qqvv_vva_mva.Histo1D(("qqvv_hvva_vv_m_mva", "", *bins_m), "vv_m"))
    hists.append(df_qqvv_vva_mva.Histo1D(("qqvv_hvva_vva_m_mva", "", *bins_m), "vva_m"))
    hists.append(df_qqvv_vva_mva.Histo1D(("qqvv_hvva_vva_p_mva", "", *bins_m), "vva_p"))
    hists.append(df_qqvv_vva_mva.Histo1D(("qqvv_hvva_vv_recoil_m_mva", "", *bins_m), "vv_recoil_m"))

    ## mass differences
    hists.append(df_qqvv_vva_mva.Histo1D(("qqvv_hvva_mass_difference_qqa_vv_mva", "", *bins_m_mass_difference), "mass_difference_qqa_vv"))
    hists.append(df_qqvv_vva_mva.Histo1D(("qqvv_hvva_mass_difference_vva_qq_mva", "", *bins_m_mass_difference), "mass_difference_vva_qq"))

    hists.append(df_qqvv_vva_mva.Histo1D(("qqvv_hvva_dr_vv_qqa_nOne", "", *bins_dr), "dr_vv_qqa"))
    hists.append(df_qqvv_vva_mva.Histo1D(("qqvv_hvva_dr_qq_vva_nOne", "", *bins_dr), "dr_qq_vva"))
    hists.append(df_qqvv_vva_mva.Histo1D(("qqvv_hvva_dr_a_qq_nOne", "", *bins_dr), "dr_a_qq"))
    hists.append(df_qqvv_vva_mva.Histo1D(("qqvv_hvva_dr_a_vv_nOne", "", *bins_dr), "dr_a_vv"))

    mva_vars = ["H_m", "H_p", "Z_m", "Z_p", "ZH_m", "ZH_p", "jet1_theta", "jet2_theta", "jet1_p", "jet2_p", "photon_p", "photon_theta", "H_recoil_m", "H_recoil_p", "Z_recoil_m", "Z_recoil_p", "mass_difference_H_Z", "dr_H_Z", "dr_a_H", "dr_a_Z", "cos_H", "cos_Z", "vv_trans", "qq_trans", "vv_long", "qq_long"]
    #return hists, weightsum, df_qqvv_vva_mva, mva_vars

    df_qqvv_vva_mva = tmva_helper_qqvv_vva_splitting_mva.run_inference(df_qqvv_vva_mva, col_name="hvva_mva_score")
    hists.append(df_qqvv_vva_mva.Histo1D(("qqvv_hvva_mva_score", "", *(1000, 0, 1)), "hvva_mva_score"))


    varbins_mva_ = [0, 0.95, 1] # 0.75 = max significance
    varbins_m_ = list(range(110, 141, 1))
    varbins_mva = array.array('d', varbins_mva_)
    varbins_m = array.array('d', varbins_m_)
    model_vva = ROOT.RDF.TH2DModel("qqvv_hvva_final_2D_mva", "", len(varbins_m_)-1, varbins_m, len(varbins_mva_)-1, varbins_mva)
    hists.append(df_qqvv_vva_mva.Histo2D(model_vva, "qq_recoil_m", "hvva_mva_score"))



    return hists, weightsum, df, []

















    ############################################################################
    ## qqqq analysis #####
    ############################################################################






    ####
    ## CUT 7: missing mass
    ####


    hists.append(df_qqqq.Histo1D(("qqqq_missingEnergy_nOne", "", *bins_m), "missingEnergy"))


    df_qqqq = df_qqqq.Define("RP_px", "FCCAnalyses::ReconstructedParticle::get_px(rps_no_photon)")
    df_qqqq = df_qqqq.Define("RP_py", "FCCAnalyses::ReconstructedParticle::get_py(rps_no_photon)")
    df_qqqq = df_qqqq.Define("RP_pz","FCCAnalyses::ReconstructedParticle::get_pz(rps_no_photon)")
    df_qqqq = df_qqqq.Define("RP_e", "FCCAnalyses::ReconstructedParticle::get_e(rps_no_photon)")
    df_qqqq = df_qqqq.Define("RP_m", "FCCAnalyses::ReconstructedParticle::get_mass(rps_no_photon)")
    df_qqqq = df_qqqq.Define("RP_q", "FCCAnalyses::ReconstructedParticle::get_charge(rps_no_photon)")
    df_qqqq = df_qqqq.Define("pseudo_jets", "FCCAnalyses::JetClusteringUtils::set_pseudoJets(RP_px, RP_py, RP_pz, RP_e)")

    df_qqqq = df_qqqq.Define("clustered_jets", "JetClustering::clustering_ee_kt(2, 4, 0, 10)(pseudo_jets)")
    df_qqqq = df_qqqq.Define("jets", "FCCAnalyses::JetClusteringUtils::get_pseudoJets(clustered_jets)")
    df_qqqq = df_qqqq.Define("jetconstituents", "FCCAnalyses::JetClusteringUtils::get_constituents(clustered_jets)")
    df_qqqq = df_qqqq.Define("jets_n", "jetconstituents.size()")
    df_qqqq = df_qqqq.Define("jets_e", "FCCAnalyses::JetClusteringUtils::get_e(jets)")
    df_qqqq = df_qqqq.Define("jets_px", "FCCAnalyses::JetClusteringUtils::get_px(jets)")
    df_qqqq = df_qqqq.Define("jets_py", "FCCAnalyses::JetClusteringUtils::get_py(jets)")
    df_qqqq = df_qqqq.Define("jets_pz", "FCCAnalyses::JetClusteringUtils::get_pz(jets)")
    df_qqqq = df_qqqq.Define("jets_m", "FCCAnalyses::JetClusteringUtils::get_m(jets)")


    df_qqqq = df_qqqq.Define("dmerge_01", "FCCAnalyses::JetClusteringUtils::get_exclusive_dmerge(clustered_jets, 0)")
    df_qqqq = df_qqqq.Define("dmerge_12", "FCCAnalyses::JetClusteringUtils::get_exclusive_dmerge(clustered_jets, 1)")
    df_qqqq = df_qqqq.Define("dmerge_23", "FCCAnalyses::JetClusteringUtils::get_exclusive_dmerge(clustered_jets, 2)")
    df_qqqq = df_qqqq.Define("dmerge_34", "FCCAnalyses::JetClusteringUtils::get_exclusive_dmerge(clustered_jets, 3)")
    df_qqqq = df_qqqq.Define("dmerge_45", "FCCAnalyses::JetClusteringUtils::get_exclusive_dmerge(clustered_jets, 4)")
    hists.append(df_qqqq.Histo1D(("qqqq_dmerge_01", "", *bins_merge), "dmerge_01"))
    hists.append(df_qqqq.Histo1D(("qqqq_dmerge_12", "", *bins_merge), "dmerge_12"))
    hists.append(df_qqqq.Histo1D(("qqqq_dmerge_23", "", *bins_merge), "dmerge_23"))
    hists.append(df_qqqq.Histo1D(("qqqq_dmerge_34", "", *bins_merge), "dmerge_34"))
    hists.append(df_qqqq.Histo1D(("qqqq_dmerge_45", "", *bins_merge), "dmerge_45"))
    

    hists.append(df_qqqq.Histo1D(("qqqq_dmerge_34_nOne", "", *bins_merge), "dmerge_34"))
    #df_qqqq = df_qqqq.Filter("dmerge_34 > 575")
    df_qqqq = df_qqqq.Filter("dmerge_34 > 800")
    hists.append(df_qqqq.Histo1D(("qqqq_cutFlow", "", *bins_count), "cut8"))

    hists.append(df_qqqq.Histo1D(("qqqq_dmerge_45_nOne", "", *bins_merge), "dmerge_45"))
    #df_qqqq = df_qqqq.Filter("dmerge_45 > 0")
    #hists.append(df_qqqq.Histo1D(("qqqq_cutFlow", "", *bins_count), "cut10"))



    #### 
    ## CUT8: JET CLUSTERING
    ####
    # require good jet clustering
    hists.append(df_qqqq.Histo1D(("qqqq_jets_n_nOne", "", *bins_count), "jets_n"))
    df_qqqq = df_qqqq.Filter("jets_n == 4")
    hists.append(df_qqqq.Histo1D(("qqqq_cutFlow", "", *bins_count), "cut9"))

    ### cuts on jet momenta/constituents?
    df_qqqq = df_qqqq.Define("jets_tlv", "FCCAnalyses::makeLorentzVectors(jets_px, jets_py, jets_pz, jets_e)")
    df_qqqq = df_qqqq.Define("jet1_p", "jets_tlv[0].P()")
    df_qqqq = df_qqqq.Define("jet2_p", "jets_tlv[1].P()")
    df_qqqq = df_qqqq.Define("jet3_p", "jets_tlv[2].P()")
    df_qqqq = df_qqqq.Define("jet4_p", "jets_tlv[3].P()")
    df_qqqq = df_qqqq.Define("jet1_theta", "jets_tlv[0].Theta()")
    df_qqqq = df_qqqq.Define("jet2_theta", "jets_tlv[1].Theta()")
    df_qqqq = df_qqqq.Define("jet3_theta", "jets_tlv[2].Theta()")
    df_qqqq = df_qqqq.Define("jet4_theta", "jets_tlv[3].Theta()")
    df_qqqq = df_qqqq.Define("jet1_nc", "jetconstituents[0].size()")
    df_qqqq = df_qqqq.Define("jet2_nc", "jetconstituents[1].size()")
    df_qqqq = df_qqqq.Define("jet3_nc", "jetconstituents[2].size()")
    df_qqqq = df_qqqq.Define("jet4_nc", "jetconstituents[3].size()")




    ####
    ## CUT9: JET ID
    ####
    # jet constituents
    hists.append(df_qqqq.Histo1D(("qqqq_jet1_nc_nOne", "", *(100, 0, 100)), "jet1_nc"))
    hists.append(df_qqqq.Histo1D(("qqqq_jet2_nc_nOne", "", *(100, 0, 100)), "jet2_nc"))
    hists.append(df_qqqq.Histo1D(("qqqq_jet3_nc_nOne", "", *(100, 0, 100)), "jet3_nc"))
    hists.append(df_qqqq.Histo1D(("qqqq_jet4_nc_nOne", "", *(100, 0, 100)), "jet4_nc"))



    df_qqqq = df_qqqq.Filter("jet1_nc > 5 && jet2_nc > 5 && jet3_nc > 5 && jet4_nc > 5")

    hists.append(df_qqqq.Histo1D(("qqqq_jet1_p_nOne", "", *bins_m), "jet1_p"))
    hists.append(df_qqqq.Histo1D(("qqqq_jet2_p_nOne", "", *bins_m), "jet2_p"))
    hists.append(df_qqqq.Histo1D(("qqqq_jet3_p_nOne", "", *bins_m), "jet3_p"))
    hists.append(df_qqqq.Histo1D(("qqqq_jet4_p_nOne", "", *bins_m), "jet4_p"))

    #df_qqqq = df_qqqq.Filter("jet1_p < 75 && jet2_p < 75 && jet3_p > 30 && jet4_p > 30")
    df_qqqq = df_qqqq.Filter("jet1_p > 5 && jet2_p > 5 && jet3_p > 5 && jet4_p > 5")
    hists.append(df_qqqq.Histo1D(("qqqq_cutFlow", "", *bins_count), "cut10"))

    # invariant mass of 4 jet system
    df_qqqq = df_qqqq.Define("jets_tot", "jets_tlv[0] + jets_tlv[1] + jets_tlv[2] + jets_tlv[3]")
    df_qqqq = df_qqqq.Define("jets_tot_m", "jets_tot.M()")
    df_qqqq = df_qqqq.Define("jets_tot_p", "jets_tot.P()")
    
    hists.append(df_qqqq.Histo1D(("qqqq_jets_tot_m_nOne", "", *bins_m), "jets_tot_m"))
    hists.append(df_qqqq.Histo1D(("qqqq_jets_tot_p_nOne", "", *bins_m), "jets_tot_p"))

    # pair jets to 2 Z bosons
    df_qqqq = df_qqqq.Define("z_pairs", "FCCAnalyses::pairing_Z(jets_tlv)")
    df_qqqq = df_qqqq.Define("Z1", "z_pairs[0]")
    df_qqqq = df_qqqq.Define("Z2", "z_pairs[1]")

    df_qqqq = df_qqqq.Define("Z1_m", "Z1.M()")
    df_qqqq = df_qqqq.Define("Z2_m", "Z2.M()")
    df_qqqq = df_qqqq.Define("Z1_p", "Z1.P()")
    df_qqqq = df_qqqq.Define("Z2_p", "Z2.P()")
    
    hists.append(df_qqqq.Histo1D(("qqqq_Z1_m", "", *bins_m), "Z1_m"))
    hists.append(df_qqqq.Histo1D(("qqqq_Z2_m", "", *bins_m), "Z2_m"))
    hists.append(df_qqqq.Histo1D(("qqqq_Z1_p", "", *bins_m), "Z1_p"))
    hists.append(df_qqqq.Histo1D(("qqqq_Z2_p", "", *bins_m), "Z2_p"))


    # pair the photon to the Z boson s
    df_qqqq = df_qqqq.Define("h_z_pairs", "FCCAnalyses::pairing_H_Z_a(ecm, Z1, Z2, photon_tlv)")
    df_qqqq = df_qqqq.Define("H", "h_z_pairs[0]") # Higgs
    df_qqqq = df_qqqq.Define("Z", "h_z_pairs[1]") # associated Z
    df_qqqq = df_qqqq.Define("ZH", "h_z_pairs[2]") # Z from Higgs

    df_qqqq = df_qqqq.Define("H_m", "H.M()")
    df_qqqq = df_qqqq.Define("H_p", "H.P()")
    df_qqqq = df_qqqq.Define("Z_m", "Z.M()")
    df_qqqq = df_qqqq.Define("Z_p", "Z.P()")
    df_qqqq = df_qqqq.Define("ZH_m", "ZH.M()")
    df_qqqq = df_qqqq.Define("ZH_p", "ZH.P()")


    hists.append(df_qqqq.Histo1D(("qqqq_H_m", "", *bins_m), "H_m"))
    hists.append(df_qqqq.Histo1D(("qqqq_H_p", "", *bins_m), "H_p"))
    hists.append(df_qqqq.Histo1D(("qqqq_Z_m", "", *bins_m), "Z_m"))
    hists.append(df_qqqq.Histo1D(("qqqq_Z_p", "", *bins_m), "Z_p"))
    hists.append(df_qqqq.Histo1D(("qqqq_ZH_m", "", *bins_m), "ZH_m"))
    hists.append(df_qqqq.Histo1D(("qqqq_ZH_p", "", *bins_m), "ZH_p"))


    # recoil of the Higgs
    df_qqqq = df_qqqq.Define("H_recoil", "init_tlv-H")
    df_qqqq = df_qqqq.Define("H_recoil_m", "(float)H_recoil.M()") # peaks at mZ for vvH, also for qqH (symmetric)
    df_qqqq = df_qqqq.Define("H_recoil_p", "(float)H_recoil.P()")
    
    df_qqqq = df_qqqq.Define("Z_recoil", "init_tlv-Z")
    df_qqqq = df_qqqq.Define("Z_recoil_m", "(float)Z_recoil.M()") # peaks at mZ for vvH, also for qqH (symmetric)
    df_qqqq = df_qqqq.Define("Z_recoil_p", "(float)Z_recoil.P()")
    
    hists.append(df_qqqq.Histo1D(("qqqq_H_recoil_m", "", *bins_m), "H_recoil_m"))
    hists.append(df_qqqq.Histo1D(("qqqq_H_recoil_p", "", *bins_m), "H_recoil_p"))
    hists.append(df_qqqq.Histo1D(("qqqq_Z_recoil_m", "", *bins_m), "Z_recoil_m"))
    hists.append(df_qqqq.Histo1D(("qqqq_Z_recoil_p", "", *bins_m), "Z_recoil_p"))

    df_qqqq = df_qqqq.Filter("H_recoil_m < 140 && H_recoil_m > 80")
    

    # mass differences
    df_qqqq = df_qqqq.Define("mass_difference", "H_m - Z_m")
    df_qqqq = df_qqqq.Define("mass_difference_Z1a_Z2", "(Z1+photon_tlv).M() - Z2.M()")
    df_qqqq = df_qqqq.Define("mass_difference_Z2a_Z1", "(Z2+photon_tlv).M() - Z1.M()")

    hists.append(df_qqqq.Histo1D(("qqqq_mass_difference", "", *bins_m_mass_difference), "mass_difference"))
    hists.append(df_qqqq.Histo1D(("qqqq_mass_difference_Z1a_Z2", "", *bins_m_mass_difference), "mass_difference_Z1a_Z2"))
    hists.append(df_qqqq.Histo1D(("qqqq_mass_difference_Z2a_Z1", "", *bins_m_mass_difference), "mass_difference_Z2a_Z1"))
    hists.append(df_qqqq.Histo1D(("qqqq_H_m_final", "", *bins_m), "H_m"))
    hists.append(df_qqqq.Histo1D(("qqqq_Z_recoil_m_final", "", *bins_m), "Z_recoil_m"))



    ### chi2 method
    df_qqqq = df_qqqq.Define("h_z_pairs_chi2", "FCCAnalyses::pairing_H_Z_a_chi2(ecm, Z1, Z2, photon_tlv)")
    hists.append(df_qqqq.Histo1D(("qqqq_chi2", "", *(2400, -200, 200)), "h_z_pairs_chi2"))

    df_qqqq_hz1 = df_qqqq.Filter("h_z_pairs_chi2 < 0") # H decays to Z1+a
    hists.append(df_qqqq_hz1.Histo1D(("qqqq_cutFlow", "", *bins_count), "cut11"))
    df_qqqq_hz1 = df_qqqq_hz1.Define("H_hz1", "Z1+photon_tlv") # Higgs
    df_qqqq_hz1 = df_qqqq_hz1.Alias("Z_hz1", "Z2") # associated Z
    df_qqqq_hz1 = df_qqqq_hz1.Alias("ZH_hz1", "Z1") # Z from Higgs
    
    df_qqqq_hz1 = df_qqqq_hz1.Define("H_m_hz1", "H_hz1.M()")
    df_qqqq_hz1 = df_qqqq_hz1.Define("H_p_hz1", "H_hz1.P()")
    df_qqqq_hz1 = df_qqqq_hz1.Define("Z_m_hz1", "Z_hz1.M()")
    df_qqqq_hz1 = df_qqqq_hz1.Define("Z_p_hz1", "Z_hz1.P()")
    df_qqqq_hz1 = df_qqqq_hz1.Define("ZH_m_hz1", "ZH_hz1.M()")
    df_qqqq_hz1 = df_qqqq_hz1.Define("ZH_p_hz1", "ZH_hz1.P()")

    hists.append(df_qqqq_hz1.Histo1D(("qqqq_H_m_hz1", "", *bins_m), "H_m_hz1"))
    hists.append(df_qqqq_hz1.Histo1D(("qqqq_H_p_hz1", "", *bins_m), "H_p_hz1"))
    hists.append(df_qqqq_hz1.Histo1D(("qqqq_Z_m_hz1", "", *bins_m), "Z_m_hz1"))
    hists.append(df_qqqq_hz1.Histo1D(("qqqq_Z_p_hz1", "", *bins_m), "Z_p_hz1"))
    hists.append(df_qqqq_hz1.Histo1D(("qqqq_ZH_m_hz1", "", *bins_m), "ZH_m_hz1"))
    hists.append(df_qqqq_hz1.Histo1D(("qqqq_ZH_p_hz1", "", *bins_m), "ZH_p_hz1"))


    # recoil of the Higgs
    df_qqqq_hz1 = df_qqqq_hz1.Define("H_recoil_hz1", "init_tlv-H_hz1")
    df_qqqq_hz1 = df_qqqq_hz1.Define("H_recoil_m_hz1", "(float)H_recoil_hz1.M()") # peaks at mZ for vvH, also for qqH (symmetric)
    df_qqqq_hz1 = df_qqqq_hz1.Define("H_recoil_p_hz1", "(float)H_recoil_hz1.P()")
    
    df_qqqq_hz1 = df_qqqq_hz1.Define("Z_recoil_hz1", "init_tlv-Z_hz1")
    df_qqqq_hz1 = df_qqqq_hz1.Define("Z_recoil_m_hz1", "(float)Z_recoil_hz1.M()") # peaks at mZ for vvH, also for qqH (symmetric)
    df_qqqq_hz1 = df_qqqq_hz1.Define("Z_recoil_p_hz1", "(float)Z_recoil_hz1.P()")
    
    hists.append(df_qqqq_hz1.Histo1D(("qqqq_H_recoil_m_hz1", "", *bins_m), "H_recoil_m_hz1"))
    hists.append(df_qqqq_hz1.Histo1D(("qqqq_H_recoil_p_hz1", "", *bins_m), "H_recoil_p_hz1"))
    hists.append(df_qqqq_hz1.Histo1D(("qqqq_Z_recoil_m_hz1", "", *bins_m), "Z_recoil_m_hz1"))
    hists.append(df_qqqq_hz1.Histo1D(("qqqq_Z_recoil_p_hz1", "", *bins_m), "Z_recoil_p_hz1"))


    
    
    df_qqqq_hz2 = df_qqqq.Filter("h_z_pairs_chi2 > 0") # H decays to Z2+a
    hists.append(df_qqqq_hz2.Histo1D(("qqqq_cutFlow", "", *bins_count), "cut12"))
    df_qqqq_hz2 = df_qqqq_hz2.Define("H_hz2", "Z2+photon_tlv") # Higgs
    df_qqqq_hz2 = df_qqqq_hz2.Alias("Z_hz2", "Z1") # associated Z
    df_qqqq_hz2 = df_qqqq_hz2.Alias("ZH_hz2", "Z2") # Z from Higgs
    
    df_qqqq_hz2 = df_qqqq_hz2.Define("H_m_hz2", "H_hz2.M()")
    df_qqqq_hz2 = df_qqqq_hz2.Define("H_p_hz2", "H_hz2.P()")
    df_qqqq_hz2 = df_qqqq_hz2.Define("Z_m_hz2", "Z_hz2.M()")
    df_qqqq_hz2 = df_qqqq_hz2.Define("Z_p_hz2", "Z_hz2.P()")
    df_qqqq_hz2 = df_qqqq_hz2.Define("ZH_m_hz2", "ZH_hz2.M()")
    df_qqqq_hz2 = df_qqqq_hz2.Define("ZH_p_hz2", "ZH_hz2.P()")

    hists.append(df_qqqq_hz2.Histo1D(("qqqq_H_m_hz2", "", *bins_m), "H_m_hz2"))
    hists.append(df_qqqq_hz2.Histo1D(("qqqq_H_p_hz2", "", *bins_m), "H_p_hz2"))
    hists.append(df_qqqq_hz2.Histo1D(("qqqq_Z_m_hz2", "", *bins_m), "Z_m_hz2"))
    hists.append(df_qqqq_hz2.Histo1D(("qqqq_Z_p_hz2", "", *bins_m), "Z_p_hz2"))
    hists.append(df_qqqq_hz2.Histo1D(("qqqq_ZH_m_hz2", "", *bins_m), "ZH_m_hz2"))
    hists.append(df_qqqq_hz2.Histo1D(("qqqq_ZH_p_hz2", "", *bins_m), "ZH_p_hz2"))


    # recoil of the Higgs
    df_qqqq_hz2 = df_qqqq_hz2.Define("H_recoil_hz2", "init_tlv-H_hz2")
    df_qqqq_hz2 = df_qqqq_hz2.Define("H_recoil_m_hz2", "(float)H_recoil_hz2.M()") # peaks at mZ for vvH, also for qqH (symmetric)
    df_qqqq_hz2 = df_qqqq_hz2.Define("H_recoil_p_hz2", "(float)H_recoil_hz2.P()")
    
    df_qqqq_hz2 = df_qqqq_hz2.Define("Z_recoil_hz2", "init_tlv-Z_hz2")
    df_qqqq_hz2 = df_qqqq_hz2.Define("Z_recoil_m_hz2", "(float)Z_recoil_hz2.M()") # peaks at mZ for vvH, also for qqH (symmetric)
    df_qqqq_hz2 = df_qqqq_hz2.Define("Z_recoil_p_hz2", "(float)Z_recoil_hz2.P()")
    
    hists.append(df_qqqq_hz2.Histo1D(("qqqq_H_recoil_m_hz2", "", *bins_m), "H_recoil_m_hz2"))
    hists.append(df_qqqq_hz2.Histo1D(("qqqq_H_recoil_p_hz2", "", *bins_m), "H_recoil_p_hz2"))
    hists.append(df_qqqq_hz2.Histo1D(("qqqq_Z_recoil_m_hz2", "", *bins_m), "Z_recoil_m_hz2"))
    hists.append(df_qqqq_hz2.Histo1D(("qqqq_Z_recoil_p_hz2", "", *bins_m), "Z_recoil_p_hz2"))

    return hists, weightsum, df_qqqq, []

    #df_qqqq = df_qqqq.Define("qq_costheta", "std::cos(qq.Theta())")




    #### 
    ## CUT10: qq mass
    ####
    # first select good Z candidate based on mass
    hists.append(df_qqqq.Histo1D(("qq_m_nOne", "", *bins_m), "qq_m"))
    #df_qqqq = df_qqqq.Filter("qq_m < (91.2+10) && qq_m > (91.2-10)")
    df_qqqq = df_qqqq.Filter("qq_m < 105 && qq_m > 80.0") # mass candidate
    hists.append(df_qqqq.Histo1D(("cutFlow", "", *bins_count), "cut10"))


    #### 
    ## CUT11: qq costheta
    ####
    hists.append(df_qqqq.Histo1D(("qq_costheta_nOne", "", *(100, 0 , 1)), "qq_costheta"))
    #df_qqqq = df_qqqq.Filter("qq_costheta < 0.7") # quite tight cut
    #df_qqqq = df_qqqq.Filter("qq_costheta < 0.8") # quite tight cut
    #hists.append(df_qqqq.Histo1D(("cutFlow", "", *bins_count), "cut11"))



    # now we need to know whether the qq comes from the Higgs Z (hqq) or from the associated Z (zqq)


    # qqa system
    df_qqqq = df_qqqq.Define("qqa", "qq_tlv+photon_tlv")
    df_qqqq = df_qqqq.Define("qqa_m", "(float)qqa.M()") # peaks at mH for vvH (sharp peak), broader peak at mH for qqH
    df_qqqq = df_qqqq.Define("qqa_p", "(float)qqa.P()")

    # recoil of qqa system
    df_qqqq = df_qqqq.Define("qqa_recoil", "init_tlv-qqa")
    df_qqqq = df_qqqq.Define("qqa_recoil_m", "(float)qqa_recoil.M()") # peaks at mZ for vvH, also for qqH (symmetric)
    df_qqqq = df_qqqq.Define("qqa_recoil_p", "(float)qqa_recoil.P()")

    # recoil of the qq system
    df_qqqq = df_qqqq.Define("qq_recoil", "init_tlv-qq_tlv") # peaks at mH for qqH, what for vvH?
    df_qqqq = df_qqqq.Define("qq_recoil_m", "(float)qq_recoil.M()")
    
    
    # vva system
    df_qqqq = df_qqqq.Define("missingMass_vec", "ROOT::Math::PxPyPzMVector(missingEnergy_tlv.Px(), missingEnergy_tlv.Py(), missingEnergy_tlv.Pz(), missingMass)")
    df_qqqq = df_qqqq.Define("missingMass_tlv", "TLorentzVector ret; ret.SetPxPyPzE(missingMass_vec.Px(), missingMass_vec.Py(), missingMass_vec.Pz(), missingMass_vec.E()); return ret;")
    df_qqqq = df_qqqq.Define("vva", "missingMass_tlv+photon_tlv")
    df_qqqq = df_qqqq.Define("vva_m", "(float)vva.M()")
    df_qqqq = df_qqqq.Define("vva_p", "(float)vva.P()")

    # recoil of vva system
    df_qqqq = df_qqqq.Define("vva_recoil", "init_tlv-vva")
    df_qqqq = df_qqqq.Define("vva_recoil_m", "vva_recoil.M()")
    df_qqqq = df_qqqq.Define("vva_recoil_p", "vva_recoil.P()")

    # recoil of vv system
    df_qqqq = df_qqqq.Define("vv_m", "(float)missingMass_tlv.M()")
    df_qqqq = df_qqqq.Define("vv_recoil", "init_tlv-missingMass_tlv") # higgs
    df_qqqq = df_qqqq.Define("vv_recoil_m", "(float)vv_recoil.M()")


    #### 
    ## CUT11: cut on recoil of the qqa system (qqa_recoil_m peaks the same for vvqqa or qqvva --> generic cut)
    ####
    hists.append(df_qqqq.Histo1D(("qqa_recoil_m_nOne", "", *bins_m), "qqa_recoil_m"))
    df_qqqq = df_qqqq.Filter("qqa_recoil_m > 80 && qqa_recoil_m < 100") ## significance
    hists.append(df_qqqq.Histo1D(("cutFlow", "", *bins_count), "cut11"))

    hists.append(df_qqqq.Histo1D(("vva_recoil_m_nOne", "", *bins_m), "vva_recoil_m"))


    # now we need to split the events based on vvqqa or qqvva --> mass separation
    hists.append(df_qqqq.Histo1D(("qq_p", "", *bins_m), "qq_p"))
    hists.append(df_qqqq.Histo1D(("qqa_m", "", *bins_m), "qqa_m"))
    hists.append(df_qqqq.Histo1D(("qqa_p", "", *bins_m), "qqa_p"))
    hists.append(df_qqqq.Histo1D(("qq_recoil_m", "", *bins_m), "qq_recoil_m"))

    hists.append(df_qqqq.Histo1D(("vva_m", "", *bins_m), "vva_m"))
    hists.append(df_qqqq.Histo1D(("vva_p", "", *bins_m), "vva_p"))
    hists.append(df_qqqq.Histo1D(("vv_recoil_m", "", *bins_m), "vv_recoil_m"))

    # angular variables
    
    # angle between Z and H
    df_qqqq = df_qqqq.Define("dr_vv_qqa", "missingMass_tlv.DeltaR(qqa)")
    df_qqqq = df_qqqq.Define("dr_qq_vva", "qq_tlv.DeltaR(vva)")
    
    # angle between photon and vv, qq
    df_qqqq = df_qqqq.Define("dr_a_qq", "photon_tlv.DeltaR(qq_tlv)")
    df_qqqq = df_qqqq.Define("dr_a_vv", "photon_tlv.DeltaR(missingMass_tlv)")

    # mass differences
    df_qqqq = df_qqqq.Define("mass_difference_qqa_vv", "qqa_m - vv_m")
    df_qqqq = df_qqqq.Define("mass_difference_vva_qq", "vva_m - qq_m")

    # cosine angles
    df_qqqq = df_qqqq.Define("cos_qqa", "abs(cos(qqa.Theta()))")
    df_qqqq = df_qqqq.Define("cos_qq", "abs(cos(qq_tlv.Theta()))")
    df_qqqq = df_qqqq.Define("cos_vva", "abs(cos(vva.Theta()))")

    #df_qqqq = df_qqqq.Define("chi2", "std::pow(za_m-125, 2) + 0.0*std::pow(za_p-50, 2) - 0.0*std::pow(qq_p-50, 2) - std::pow(recoil_z_m-125, 2)")
    df_qqqq = df_qqqq.Define("chi2", "std::pow(qqa_m-125, 2) + 0.0*std::pow(vv_recoil_m-125, 2) - std::pow(qq_recoil_m-125, 2) - std::pow(vva_m-125, 2)")
    # adding qq_p and za_p makes it worse
    hists.append(df_qqqq.Histo1D(("chi2", "", *(2400, -200, 200)), "chi2"))

























    return hists, weightsum, df_qqvv, []



if treemaker:
    mva_variables = []
    class RDFanalysis():
        def analysers(df):
            global mva_variables
            hists, weightsum, df, mva_vars = build_graph_za(df, "")
            mva_variables = mva_vars
            return df

        # define output branches to be saved
        def output():
            #branchList = ["qqa_m", "qq_m", "vv_recoil_m", "qq_recoil_m", "vva_m", "qqa_p", "qq_p"]
            #branchList = ["qqa_m", "qq_m", "vv_recoil_m", "qq_recoil_m", "vva_m", "qqa_p", "qq_p", "qqa_recoil_m", "vv_m", "vva_p", "dr_vv_qqa", "dr_qq_vva", "dr_a_qq" ,"dr_a_vv", "mass_difference_qqa_vv", "mass_difference_vva_qq", "cos_qqa", "cos_qq", "cos_vva", "photon_p", "photon_theta"]
 
            
            
            # qqqqa
            #branchList = vars = ["H_m", "H_p", "Z_m", "Z_p", "ZH_m", "ZH_p", "jet1_theta", "jet2_theta", "jet3_theta", "jet4_theta", "jet1_p", "jet2_p", "jet3_p", "jet4_p", "photon_p", "photon_theta", "H_recoil_m", "H_recoil_p", "Z_recoil_m", "Z_recoil_p", "mass_difference"]
            #branchList = ["H_m", "H_p", "Z_m", "Z_p", "ZH_m", "ZH_p", "dr_H_Z", "dr_a_H", "dr_a_Z", "mass_difference_H_Z", "cos_H", "cos_Z", "photon_p", "photon_theta", "jet1_theta", "jet2_theta", "jet1_p", "jet2_p", "H_recoil_m", "H_recoil_p", "Z_recoil_m", "Z_recoil_p"]
            
            ## qqvv chi2 splitting
            #branchList = ["H_m", "H_p", "Z_m", "Z_p", "ZH_m", "ZH_p", "jet1_theta", "jet2_theta", "jet1_p", "jet2_p", "photon_p", "photon_theta", "H_recoil_m", "H_recoil_p", "Z_recoil_m", "Z_recoil_p", "mass_difference_H_Z", "dr_H_Z", "dr_a_H", "dr_a_Z", "cos_H", "cos_Z"]
            global mva_variables
            return mva_variables

else:
    def build_graph(df, dataset):
        hists, weightsum, df, mva_vars = build_graph_za(df, dataset)
        return hists, weightsum
