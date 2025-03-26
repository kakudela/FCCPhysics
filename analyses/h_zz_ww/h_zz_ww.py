
import ROOT
import array
ROOT.TH1.SetDefaultSumw2(ROOT.kTRUE)
from addons.TMVAHelper.TMVAHelper import TMVAHelperXGB

ecm = 240 # 240 365
treemaker = False

# list of all processes
fraction = 0.01


if ecm == 240:
    processList = {
        'wzp6_ee_nunuH_HWW_ecm240':         {'fraction':1},
        'wzp6_ee_eeH_HWW_ecm240':           {'fraction':1},
        'wzp6_ee_tautauH_HWW_ecm240':       {'fraction':1},
        'wzp6_ee_mumuH_HWW_ecm240':         {'fraction':1},
        'wzp6_ee_ccH_HWW_ecm240':           {'fraction':1},
        'wzp6_ee_bbH_HWW_ecm240':           {'fraction':1},
        'wzp6_ee_qqH_HWW_ecm240':           {'fraction':1},
        'wzp6_ee_ssH_HWW_ecm240':           {'fraction':1},

        'wzp6_ee_nunuH_HZZ_ecm240':         {'fraction':1},
        'wzp6_ee_eeH_HZZ_ecm240':           {'fraction':1},
        'wzp6_ee_tautauH_HZZ_ecm240':       {'fraction':1},
        'wzp6_ee_mumuH_HZZ_ecm240':         {'fraction':1},
        'wzp6_ee_ccH_HZZ_ecm240':           {'fraction':1},
        'wzp6_ee_bbH_HZZ_ecm240':           {'fraction':1},
        'wzp6_ee_qqH_HZZ_ecm240':           {'fraction':1},
        'wzp6_ee_ssH_HZZ_ecm240':           {'fraction':1},

        'wz3p6_ee_uu_ecm240':               {'fraction':fraction},
        'wz3p6_ee_dd_ecm240':               {'fraction':fraction},
        'wz3p6_ee_cc_ecm240':               {'fraction':fraction},
        'wz3p6_ee_ss_ecm240':               {'fraction':fraction},
        'wz3p6_ee_bb_ecm240':               {'fraction':fraction},
        #'wz3p6_ee_gammagamma_ecm240':       {'fraction':fraction},
        #'wz3p6_ee_tautau_ecm240':           {'fraction':fraction},
        #'wz3p6_ee_mumu_ecm240':             {'fraction':fraction},
        #'wz3p6_ee_ee_Mee_30_150_ecm240':    {'fraction':fraction},
        #'wz3p6_ee_nunu_ecm240':             {'fraction':fraction},

        'p8_ee_ZZ_ecm240':             {'fraction':fraction},
        'p8_ee_WW_ecm240':             {'fraction':fraction},
    }





inputDir = "/ceph/submit/data/group/fcc/ee/generation/DelphesEvents/winter2023/IDEA/"
procDict = "/ceph/submit/data/group/fcc/ee/generation/DelphesEvents/winter2023/IDEA/samplesDict.json"

# additional/custom C++ functions
includePaths = ["../../functions/functions.h", "../../functions/functions_gen.h", "../h_zh/utils.h", "utils.h"]


# output directory
if treemaker:
    outputDir   = f"/ceph/submit/data/group/fcc/ee/analyses/h_zz_ww/treemaker/ecm{ecm}/"
else:
    outputDir   = f"output/h_zz_ww/histmaker/ecm{ecm}/"

# optional: ncpus, default is 4, -1 uses all cores available
nCPUS       = 32

# scale the histograms with the cross-section and integrated luminosity
doScale = True
intLumi = 10.8e6 if ecm == 240 else 3e6

# define histograms

# define histograms
bins_m = (400, 0, 400) # 100 MeV bins
bins_p_mu = (2000, 0, 200) # 100 MeV bins
bins_m_ll = (2000, 0, 200) # 100 MeV bins
bins_p_ll = (200, 0, 200) # 1 GeV bins
bins_recoil = (20000, 0, 200) # 10 MeV bins 
bins_recoil_fine = (20000, 120, 140) # 1 MeV bins 
bins_cosThetaMiss = (10000, 0, 1)

bins_theta = (500, 0, 5)
bins_phi = (500, -5, 5)
bins_aco = (1000, -10, 10)

bins_count = (50, 0, 50)
bins_pdgid = (60, -30, 30)
bins_charge = (10, -5, 5)
bins_iso = (500, 0, 5)
bins_dR = (1000, 0, 10)

bins_cat = (10, 0, 10)
bins_resolution = (10000, 0.95, 1.05)

bins_m_fine = (100, 120, 130) # 100 MeV bins
bins_cos_abs = (100, 0, 1)

bins_merge = (50000, 0, 50000)
bins_merge_sq = (500, 0, 500)


ROOT.EnableImplicitMT(nCPUS) # hack to deal correctly with TMVAHelperXGB  # bdt_model_new_0p1_WZqq bdt_model_new_0p1 bdt_model_0p1_inv bdt_model_0p1
## original: output/h_zh_hadronic/training/bdt_model_WW_Zg_thrust_reduced.root
## before: /ceph/submit/data/group/fcc/ee/analyses/zh/hadronic/training/bdt_model_final.root
#tmva_helper = TMVAHelperXGB("/ceph/submit/data/group/fcc/ee/analyses/zh/hadronic/training/bdt_model_ecm240.root", "bdt_model") # read the XGBoost training
#tmva_helper = TMVAHelperXGB("output/h_zh_hadronic/training/bdt_model_WW_Zg_thrust_reduced_withInv.root", "bdt_model") # read the XGBoost training





def exclusive_clustering(df, njets):
    if njets==0: # inclusive
        df = df.Define("clustered_jets_N0", "JetClustering::clustering_kt(0.6, 0, 5, 1, 0)(pseudo_jets)")
    else:
        df = df.Define(f"clustered_jets_N{njets}", f"JetClustering::clustering_ee_kt(2, {njets}, 0, 10)(pseudo_jets)")
        #df = df.Define(f"clustered_jets_N{njets}", f"JetClustering::clustering_ee_genkt(0.5, 2, {njets}, 0, 10, 0)(pseudo_jets)")
    df = df.Define(f"jets_N{njets}", f"FCCAnalyses::JetClusteringUtils::get_pseudoJets(clustered_jets_N{njets})")
    df = df.Define(f"njets_N{njets}", f"jets_N{njets}.size()")
    df = df.Define(f"jetconstituents_N{njets}", f"FCCAnalyses::JetClusteringUtils::get_constituents(clustered_jets_N{njets})")
    df = df.Define(f"jets_e_N{njets}", f"FCCAnalyses::JetClusteringUtils::get_e(jets_N{njets})")
    df = df.Define(f"jets_px_N{njets}", f"FCCAnalyses::JetClusteringUtils::get_px(jets_N{njets})")
    df = df.Define(f"jets_py_N{njets}", f"FCCAnalyses::JetClusteringUtils::get_py(jets_N{njets})")
    df = df.Define(f"jets_pz_N{njets}", f"FCCAnalyses::JetClusteringUtils::get_pz(jets_N{njets})")
    df = df.Define(f"jets_m_N{njets}", f"FCCAnalyses::JetClusteringUtils::get_m(jets_N{njets})")
    df = df.Define(f"jets_rp_N{njets}", f"FCCAnalyses::jets2rp(jets_px_N{njets}, jets_py_N{njets}, jets_pz_N{njets}, jets_e_N{njets}, jets_m_N{njets})")
    df = df.Define(f"jets_rp_cand_N{njets}", f"FCCAnalyses::select_jets(jets_rp_N{njets}, jetconstituents_N{njets}, {njets}, ReconstructedParticles)") # reduces potentially the jet multiplicity
    df = df.Define(f"njets_cand_N{njets}", f"jets_rp_cand_N{njets}.size()")
    return df

def define_variables(df, njets):
    #if njets == 0:
    #df = df.Filter(f"if(njets_cand_N{njets} < 2) std::cout << \"FILTER={njets} \" << njets_cand_N{njets}  <<  std::endl; return njets_cand_N{njets} >= 2") #### tricky need to fix it

    df = df.Define(f"zbuilder_result_N{njets}", f"FCCAnalyses::resonanceBuilder_mass_recoil_hadronic(91.2, 125, 0.0, ecm)(jets_rp_cand_N{njets})")
    df = df.Define(f"zqq_N{njets}", f"ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>{{zbuilder_result_N{njets}[0]}}") # the Z
    df = df.Define(f"zqq_jets_N{njets}", f"ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>{{zbuilder_result_N{njets}[1],zbuilder_result_N{njets}[2]}}") # the leptons
    df = df.Define(f"zqq_m_N{njets}", f"FCCAnalyses::ReconstructedParticle::get_mass(zqq_N{njets})[0]")
    df = df.Define(f"zqq_p_N{njets}", f"FCCAnalyses::ReconstructedParticle::get_p(zqq_N{njets})[0]")
    df = df.Define(f"zqq_recoil_N{njets}", f"FCCAnalyses::ReconstructedParticle::recoilBuilder(ecm)(zqq_N{njets})")
    df = df.Define(f"zqq_recoil_m_N{njets}", f"FCCAnalyses::ReconstructedParticle::get_mass(zqq_recoil_N{njets})[0]")


    return df







def build_graph_hwwzz(df, dataset):

    hists = []

    df = df.Define("ecm", "240" if ecm == 240 else "365")
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

    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut0")) ## all events
    if "HWW" in dataset: # remove muons/electrons from inclusive WW ## ONLY HADRONIC
        df = df.Define("ww_hadronic", "FCCAnalyses::is_ww_hadronic(Particle, Particle1)")
        #df = df.Filter("ww_hadronic")

    if "HZZ" in dataset: # remove muons/electrons from inclusive ZZ ## ONLY HADRONIC
        df = df.Define("zz_hadronic", "FCCAnalyses::is_zz_hadronic(Particle, Particle1)")
        #df = df.Filter("zz_hadronic")


    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut1")) ## hadronic events








    ### PRESELECTION
    #hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut0")) ## all events
    #hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut0")) ## all events
    #hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut0")) ## all events

    # muons
    df = df.Alias("Muon0", "Muon#0.index")
    df = df.Define("muons_all", "FCCAnalyses::ReconstructedParticle::get(Muon0, ReconstructedParticles)")
    df = df.Define("muons_all_p", "FCCAnalyses::ReconstructedParticle::get_p(muons_all)")
    df = df.Define("muons_all_theta", "FCCAnalyses::ReconstructedParticle::get_theta(muons_all)")
    df = df.Define("muons_all_phi", "FCCAnalyses::ReconstructedParticle::get_phi(muons_all)")
    df = df.Define("muons_all_q", "FCCAnalyses::ReconstructedParticle::get_charge(muons_all)")
    df = df.Define("muons_all_no", "FCCAnalyses::ReconstructedParticle::get_n(muons_all)")

    df = df.Define("muons", "FCCAnalyses::ReconstructedParticle::sel_p(10)(muons_all)")
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

    df = df.Define("electrons", "FCCAnalyses::ReconstructedParticle::sel_p(10)(electrons_all)")
    df = df.Define("electrons_p", "FCCAnalyses::ReconstructedParticle::get_p(electrons)")
    df = df.Define("electrons_theta", "FCCAnalyses::ReconstructedParticle::get_theta(electrons)")
    df = df.Define("electrons_phi", "FCCAnalyses::ReconstructedParticle::get_phi(electrons)")
    df = df.Define("electrons_q", "FCCAnalyses::ReconstructedParticle::get_charge(electrons)")
    df = df.Define("electrons_no", "FCCAnalyses::ReconstructedParticle::get_n(electrons)")

    ## VETO leptons
    df = df.Filter("muons_no == 0 && electrons_no == 0")
    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut2"))



    ## MISSING MOMENTUM
    df = df.Define("missingEnergy_rp", "FCCAnalyses::missingEnergy(ecm, ReconstructedParticles)")
    df = df.Define("missingEnergy", "missingEnergy_rp[0].energy")
    df = df.Define("missingEnergy_p", "FCCAnalyses::ReconstructedParticle::get_p(missingEnergy_rp)[0]")
    hists.append(df.Histo1D(("missingEnergy_nOne", "", *bins_m), "missingEnergy_p"))
    df = df.Filter("missingEnergy_p < 20")
    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut3"))



    ## CLUSTERING
    df = df.Alias("rps_sel", "ReconstructedParticles")
    df = df.Define("RP_px", "FCCAnalyses::ReconstructedParticle::get_px(rps_sel)")
    df = df.Define("RP_py", "FCCAnalyses::ReconstructedParticle::get_py(rps_sel)")
    df = df.Define("RP_pz","FCCAnalyses::ReconstructedParticle::get_pz(rps_sel)")
    df = df.Define("RP_e", "FCCAnalyses::ReconstructedParticle::get_e(rps_sel)")
    df = df.Define("RP_m", "FCCAnalyses::ReconstructedParticle::get_mass(rps_sel)")
    df = df.Define("RP_q", "FCCAnalyses::ReconstructedParticle::get_charge(rps_sel)")
    df = df.Define("pseudo_jets", "FCCAnalyses::JetClusteringUtils::set_pseudoJets(RP_px, RP_py, RP_pz, RP_e)")

    df = df.Define("clustered_jets", "JetClustering::clustering_ee_kt(2, 6, 0, 10)(pseudo_jets)")
    df = df.Define("jets", "FCCAnalyses::JetClusteringUtils::get_pseudoJets(clustered_jets)")
    df = df.Define("njets", "jets.size()")
    df = df.Define("jetconstituents", "FCCAnalyses::JetClusteringUtils::get_constituents(clustered_jets)")
    df = df.Define("jets_e", "FCCAnalyses::JetClusteringUtils::get_e(jets)")
    df = df.Define("jets_px", "FCCAnalyses::JetClusteringUtils::get_px(jets)")
    df = df.Define("jets_py", "FCCAnalyses::JetClusteringUtils::get_py(jets)")
    df = df.Define("jets_pz", "FCCAnalyses::JetClusteringUtils::get_pz(jets)")
    df = df.Define("jets_m_", "FCCAnalyses::JetClusteringUtils::get_m(jets)")
    #df = df.Define("jets_rp_", "FCCAnalyses::jets2rp(jets_px, jets_py, jets_pz, jets_e, jets_m)")
    #df = df.Define("njets_cand", "jets_rp_cand.size()")
    
    


    hists.append(df.Histo1D(("njets_nOne", "", *bins_count), "njets"))
    df = df.Filter("njets == 6")
    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut4"))





    ## VISIBLE ENERGY
    df = df.Define("visibleEnergy", "FCCAnalyses::visibleEnergy(ReconstructedParticles)")
    hists.append(df.Histo1D(("visibleEnergy_nOne", "", *bins_m), "visibleEnergy"))
    df = df.Filter("visibleEnergy > 115")
    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut5"))
    
    ## VISIBLE MASS
    df = df.Define("visibleMass", "FCCAnalyses::visibleMass(ReconstructedParticles)")
    hists.append(df.Histo1D(("visibleMass_nOne", "", *bins_m), "visibleMass"))
    df = df.Filter("visibleMass > 10")
    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut6"))




    ## cuts on dmerge
    df = df.Define("dmerge_01", "FCCAnalyses::JetClusteringUtils::get_exclusive_dmerge(clustered_jets, 0)")
    df = df.Define("dmerge_12", "FCCAnalyses::JetClusteringUtils::get_exclusive_dmerge(clustered_jets, 1)")
    df = df.Define("dmerge_23", "FCCAnalyses::JetClusteringUtils::get_exclusive_dmerge(clustered_jets, 2)")
    df = df.Define("dmerge_34", "FCCAnalyses::JetClusteringUtils::get_exclusive_dmerge(clustered_jets, 3)")
    df = df.Define("dmerge_45", "FCCAnalyses::JetClusteringUtils::get_exclusive_dmerge(clustered_jets, 4)")
    df = df.Define("dmerge_56", "FCCAnalyses::JetClusteringUtils::get_exclusive_dmerge(clustered_jets, 5)")
    df = df.Define("dmerge_67", "FCCAnalyses::JetClusteringUtils::get_exclusive_dmerge(clustered_jets, 6)")

    df = df.Define("sq_dmerge_01", "sqrt(dmerge_01)")
    df = df.Define("sq_dmerge_12", "sqrt(dmerge_12)")
    df = df.Define("sq_dmerge_23", "sqrt(dmerge_23)")
    df = df.Define("sq_dmerge_34", "sqrt(dmerge_34)")
    df = df.Define("sq_dmerge_45", "sqrt(dmerge_45)")
    df = df.Define("sq_dmerge_56", "sqrt(dmerge_56)")
    df = df.Define("sq_dmerge_67", "sqrt(dmerge_67)")

    hists.append(df.Histo1D(("dmerge_01", "", *bins_merge), "dmerge_01"))
    hists.append(df.Histo1D(("sq_dmerge_01", "", *bins_merge_sq), "sq_dmerge_01"))
    hists.append(df.Histo1D(("dmerge_12", "", *bins_merge), "dmerge_12"))
    hists.append(df.Histo1D(("sq_dmerge_12", "", *bins_merge_sq), "sq_dmerge_12"))
    hists.append(df.Histo1D(("dmerge_23", "", *bins_merge), "dmerge_23"))
    hists.append(df.Histo1D(("sq_dmerge_23", "", *bins_merge_sq), "sq_dmerge_23"))

    df = df.Filter("sq_dmerge_23 > 50")
    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut7"))


    hists.append(df.Histo1D(("dmerge_34", "", *bins_merge), "dmerge_34"))
    hists.append(df.Histo1D(("sq_dmerge_34", "", *bins_merge_sq), "sq_dmerge_34"))

    df = df.Filter("sq_dmerge_34 > 25")
    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut8"))

    hists.append(df.Histo1D(("dmerge_45", "", *bins_merge), "dmerge_45"))
    hists.append(df.Histo1D(("sq_dmerge_45", "", *bins_merge_sq), "sq_dmerge_45"))

    df = df.Filter("sq_dmerge_45 > 15")
    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut9"))

    hists.append(df.Histo1D(("dmerge_56", "", *bins_merge), "dmerge_56"))
    hists.append(df.Histo1D(("sq_dmerge_56", "", *bins_merge_sq), "sq_dmerge_56"))

    df = df.Filter("sq_dmerge_56 > 10")
    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut10"))

    hists.append(df.Histo1D(("dmerge_67", "", *bins_merge), "dmerge_67"))
    hists.append(df.Histo1D(("sq_dmerge_67", "", *bins_merge_sq), "sq_dmerge_67"))

    
    
    
    
    
    
    
    df = df.Define("init_tlv", "TLorentzVector ret; ret.SetPxPyPzE(0, 0, 0, ecm); return ret;")
    df = df.Define("jets_tlv", "FCCAnalyses::makeLorentzVectors(jets_px, jets_py, jets_pz, jets_e)")
    df = df.Define("paired_jets_WW", "FCCAnalyses::pairing_WW(jets_tlv, 80.385)")
    df = df.Define("paired_jets_ZZ", "FCCAnalyses::pairing_WW(jets_tlv, 91.19)")

    df = df.Define("W1_WW", "paired_jets_WW[0]")
    df = df.Define("W2_WW", "paired_jets_WW[1]")
    df = df.Define("Z_WW", "paired_jets_WW[2]")
    df = df.Define("H_WW", "paired_jets_WW[3]")

    df = df.Define("Z1_ZZ", "paired_jets_ZZ[0]")
    df = df.Define("Z2_ZZ", "paired_jets_ZZ[1]")
    df = df.Define("Z_ZZ", "paired_jets_ZZ[2]")
    df = df.Define("H_ZZ", "paired_jets_ZZ[3]")

    df = df.Define("W1_WW_m", "W1_WW.M()")
    df = df.Define("W2_WW_m", "W2_WW.M()")
    df = df.Define("Z_WW_m", "Z_WW.M()")
    df = df.Define("Z_WW_p", "Z_WW.P()")
    df = df.Define("H_WW_m", "H_WW.M()")
    df = df.Define("H_WW_p", "H_WW.P()")

    hists.append(df.Histo1D(("W1_WW_m", "", *bins_m), "W1_WW_m"))
    hists.append(df.Histo1D(("W2_WW_m", "", *bins_m), "W2_WW_m"))
    hists.append(df.Histo1D(("Z_WW_m", "", *bins_m), "Z_WW_m"))
    hists.append(df.Histo1D(("Z_WW_p", "", *bins_m), "Z_WW_p"))
    hists.append(df.Histo1D(("H_WW_m", "", *bins_m), "H_WW_m"))
    hists.append(df.Histo1D(("H_WW_p", "", *bins_m), "H_WW_p"))

    df = df.Define("Z1_ZZ_m", "Z1_ZZ.M()")
    df = df.Define("Z2_ZZ_m", "Z2_ZZ.M()")
    df = df.Define("Z_ZZ_m", "Z_ZZ.M()")
    df = df.Define("Z_ZZ_p", "Z_ZZ.P()")
    df = df.Define("H_ZZ_m", "H_ZZ.M()")
    df = df.Define("H_ZZ_p", "H_ZZ.P()")

    hists.append(df.Histo1D(("Z1_ZZ_m", "", *bins_m), "Z1_ZZ_m"))
    hists.append(df.Histo1D(("Z2_ZZ_m", "", *bins_m), "Z2_ZZ_m"))
    hists.append(df.Histo1D(("Z_ZZ_m", "", *bins_m), "Z_ZZ_m"))
    hists.append(df.Histo1D(("Z_ZZ_p", "", *bins_m), "Z_ZZ_p"))
    hists.append(df.Histo1D(("H_ZZ_m", "", *bins_m), "H_ZZ_m"))
    hists.append(df.Histo1D(("H_ZZ_p", "", *bins_m), "H_ZZ_p"))


    return hists, weightsum, df



















    #df = df.Filter("ReconstructedParticles.size() > 30")
    
    #### remove isolated photons from clustering
    df = df.Define("photons_all_p", "FCCAnalyses::ReconstructedParticle::get_p(photons_all)")
    df = df.Define("photons", "FCCAnalyses::sel_range(40, 95, false)(photons_all, photons_all_p)")
    df = df.Define("rps_no_photons", "FCCAnalyses::ReconstructedParticle::remove(ReconstructedParticles, photons)")


    #### remove isolated muons from clustering
    # ensure orthogonality with leptonic channel
    df = df.Define("muons_all_p", "FCCAnalyses::ReconstructedParticle::get_p(muons_all)")
    df = df.Define("muons", "FCCAnalyses::sel_range(40, 95, false)(muons_all, muons_all_p)")

    #### remove isolated electrons from clustering
    # ensure orthogonality with leptonic channel
    df = df.Define("electrons_all_p", "FCCAnalyses::ReconstructedParticle::get_p(electrons_all)")
    df = df.Define("electrons", "FCCAnalyses::sel_range(40, 95, false)(electrons_all, electrons_all_p)")


    ### remove isolated electrons from clustering
    # ensure orthogonality with leptonic channel
    df = df.Define("rps_no_photons_muons", "FCCAnalyses::ReconstructedParticle::remove(rps_no_photons, muons)")
    df = df.Define("rps_no_photons_muons_electrons", "FCCAnalyses::ReconstructedParticle::remove(rps_no_photons_muons, electrons)")


    # define PF candidates collection by removing the muons
    #df = df.Alias("rps_sel", "rps_no_photons_muons_electrons")
    df = df.Alias("rps_sel", "ReconstructedParticles")
    df = df.Define("RP_px", "FCCAnalyses::ReconstructedParticle::get_px(rps_sel)")
    df = df.Define("RP_py", "FCCAnalyses::ReconstructedParticle::get_py(rps_sel)")
    df = df.Define("RP_pz","FCCAnalyses::ReconstructedParticle::get_pz(rps_sel)")
    df = df.Define("RP_e", "FCCAnalyses::ReconstructedParticle::get_e(rps_sel)")
    df = df.Define("RP_m", "FCCAnalyses::ReconstructedParticle::get_mass(rps_sel)")
    df = df.Define("RP_q", "FCCAnalyses::ReconstructedParticle::get_charge(rps_sel)")
    df = df.Define("pseudo_jets", "FCCAnalyses::JetClusteringUtils::set_pseudoJets(RP_px, RP_py, RP_pz, RP_e)")


    # perform possible clusterings
    
    df = exclusive_clustering(df, 6)
    df = exclusive_clustering(df, 4)
    df = exclusive_clustering(df, 2)
    df = exclusive_clustering(df, 0)

    df = define_variables(df, 6)
    df = define_variables(df, 4)
    df = define_variables(df, 2)
    df = define_variables(df, 0)


    df = df.Define("zqq", "std::vector<Vec_rp> r = {zqq_N0, zqq_N2, zqq_N4, zqq_N6}; return r;")
    df = df.Define("zqq_jets", "std::vector<Vec_rp> r = {jets_rp_cand_N0, jets_rp_cand_N2, jets_rp_cand_N4, jets_rp_cand_N6}; return r;")
    df = df.Define("zqq_m", "Vec_f r = {zqq_m_N0, zqq_m_N2, zqq_m_N4, zqq_m_N6}; return r;")
    df = df.Define("zqq_p", "Vec_f r = {zqq_p_N0, zqq_p_N2, zqq_p_N4, zqq_p_N6}; return r;")
    df = df.Define("zqq_recoil_m", "Vec_f r = {zqq_recoil_m_N0, zqq_recoil_m_N2, zqq_recoil_m_N4, zqq_recoil_m_N6}; return r;")
    df = df.Define("njets", "Vec_i r = {(int)njets_cand_N0, (int)njets_cand_N2, (int)njets_cand_N4, (int)njets_cand_N6}; return r;")
    df = df.Define("njets_target", "Vec_i r = {0, 2, 4, 6}; return r;")
    df = df.Define("best_clustering_idx", "FCCAnalyses::best_clustering_idx(zqq_m, zqq_p, zqq_recoil_m, njets, njets_target)")

    hists.append(df.Histo1D(("best_clustering_idx_nosel", "", *(15, -5, 10)), "best_clustering_idx"))

    ## njets for inclusive clustering
    hists.append(df.Histo1D(("njets_inclusive", "", *bins_count), "njets_cand_N0"))
    df_incl = df.Filter("best_clustering_idx == 0")
    hists.append(df_incl.Histo1D(("njets_inclusive_sel", "", *bins_count), "njets_cand_N0"))




    df = df.Filter("best_clustering_idx >= 0")
    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut2")) ## after clustering

    df = df.Define("zqq_best", "zqq[best_clustering_idx]")
    df = df.Define("zqq_jets_best", "zqq_jets[best_clustering_idx]")
    df = df.Define("zqq_m_best", "zqq_m[best_clustering_idx]")
    df = df.Define("zqq_recoil_m_best", "zqq_recoil_m[best_clustering_idx]")
    df = df.Define("zqq_p_best", "zqq_p[best_clustering_idx]")
    df = df.Define("zqq_jets_p_best", "FCCAnalyses::ReconstructedParticle::get_p(zqq_jets_best)")
    df = df.Define("zqq_jets_theta_best", "FCCAnalyses::ReconstructedParticle::get_theta(zqq_jets_best)")
    df = df.Define("zqq_jets_costheta_best", "Vec_f ret; for(auto & theta: zqq_jets_theta_best) ret.push_back(std::abs(cos(theta))); return ret;")

    df = df.Define("z_theta", "FCCAnalyses::ReconstructedParticle::get_theta(zqq_best)")
    df = df.Define("z_costheta", "std::abs(cos(z_theta[0]))")
    hists.append(df.Histo1D(("z_costheta_nosel", "", *bins_cos_abs), "z_costheta"))
    

    hists.append(df.Histo1D(("zqq_m_best_nosel", "", *(200, 0, 200)), "zqq_m_best"))
    hists.append(df.Histo1D(("zqq_p_best_nosel", "", *(200, 0, 200)), "zqq_p_best"))
    hists.append(df.Histo1D(("zqq_recoil_m_best_nosel", "", *(200, 0, 200)), "zqq_recoil_m_best"))
    hists.append(df.Histo2D(("zqq_recoil_m_mqq_nosel", "", *((25, 100, 150)+(35, 50, 120))), "zqq_recoil_m_best", "zqq_m_best"))

    ## jet kinematics
    df = df.Define("leading_idx", "(zqq_jets_p_best[0] > zqq_jets_p_best[1]) ? 0 : 1")
    df = df.Define("subleading_idx", "(zqq_jets_p_best[0] > zqq_jets_p_best[1]) ? 1 : 0")
    df = df.Define("leading_jet_p", "zqq_jets_p_best[leading_idx]")
    df = df.Define("subleading_jet_p", "zqq_jets_p_best[subleading_idx]")
    df = df.Define("leading_jet_costheta", "zqq_jets_costheta_best[leading_idx]")
    df = df.Define("subleading_jet_costheta", "zqq_jets_costheta_best[subleading_idx]")
    
    hists.append(df.Histo1D(("leading_jet_p", "", *(200, 0, 200)), "leading_jet_p"))
    hists.append(df.Histo1D(("subleading_jet_p", "", *(200, 0, 200)), "subleading_jet_p"))
    hists.append(df.Histo1D(("leading_jet_costheta", "", *bins_cos_abs), "leading_jet_costheta"))
    hists.append(df.Histo1D(("subleading_jet_costheta", "", *bins_cos_abs), "subleading_jet_costheta"))

    ## attempt to reconstruct WW with 4 jets
    
    df = df.Define("pairs_WW_N4", "FCCAnalyses::pair_WW_N4(jets_rp_cand_N4)")
    df = df.Define("W1", "pairs_WW_N4[0]")
    df = df.Define("W2", "pairs_WW_N4[1]")
    df = df.Define("W1_m", "W1.M()")
    df = df.Define("W2_m", "W2.M()")
    df = df.Define("W1_p", "W1.P()")
    df = df.Define("W2_p", "W2.P()")
    df = df.Define("W1_costheta", "std::abs(W1.Theta())")
    df = df.Define("W2_costheta", "std::abs(W1.Theta())")
    df = df.Define("delta_mWW", "std::sqrt((W1_m-78)*(W1_m-78) + (W2_m-78)*(W2_m-78))")

    hists.append(df.Histo1D(("delta_mWW_nosel", "", *(1000, 0, 100)), "delta_mWW"))
    hists.append(df.Histo1D(("W1_m_nosel", "", *(200, 0, 200)), "W1_m"))
    hists.append(df.Histo1D(("W2_m_nosel", "", *(200, 0, 200)), "W2_m"))
    hists.append(df.Histo1D(("W1_p_nosel", "", *(200, 0, 200)), "W1_p"))
    hists.append(df.Histo1D(("W2_p_nosel", "", *(200, 0, 200)), "W2_p"))
    hists.append(df.Histo1D(("W1_costheta_nosel", "", *bins_cos_abs), "W1_costheta"))
    hists.append(df.Histo1D(("W2_costheta_nosel", "", *bins_cos_abs), "W2_costheta"))
    ##
    '''
    df = df.Define("jets_rp_cand_N4_prime", "FCCAnalyses::energyReconstructFourJet_rp(240, jets_rp_cand_N4)")
    df = df.Define("pairs_WW_N4_prime", "FCCAnalyses::pair_WW_N4(jets_rp_cand_N4_prime)")
    df = df.Define("W1_prime", "pairs_WW_N4_prime[0]")
    df = df.Define("W2_prime", "pairs_WW_N4_prime[1]")
    df = df.Define("W1_m_prime", "W1_prime.M()")
    df = df.Define("W2_m_prime", "W2_prime.M()")
    df = df.Define("W1_p_prime", "W1_prime.M()")
    df = df.Define("W2_p_prime", "W2_prime.M()")
    df = df.Define("delta_mWW_prime", "std::sqrt((W1_m_prime-78)*(W1_m_prime-78) + (W2_m_prime-78)*(W2_m_prime-78))")

    hists.append(df.Histo1D(("delta_mWW_nosel_prime", "", *(1000, 0, 100)), "delta_mWW_prime"))
    hists.append(df.Histo1D(("W1_m_nosel_prime", "", *(200, 0, 200)), "W1_m_prime"))
    hists.append(df.Histo1D(("W2_m_nosel_prime", "", *(200, 0, 200)), "W2_m_prime"))
    hists.append(df.Histo1D(("W1_p_nosel_prime", "", *(200, 0, 200)), "W1_p_prime"))
    hists.append(df.Histo1D(("W2_p_nosel_prime", "", *(200, 0, 200)), "W2_p_prime"))
    '''
    ##############################################################
    ## START SELECTION
    ##############################################################



    ## CUT ON m(qq)
    
    hists.append(df.Histo1D(("zqq_m_best_nOne", "", *(200, 0, 200)), "zqq_m_best"))
    ##df = df.Filter("zqq_m_best > 60 && zqq_m_best < 110") # tighter at 80, but aa/mumu problem
    df = df.Filter("zqq_m_best > 20 && zqq_m_best < 140") ## loose
    #df = df.Filter("zqq_m_best > 40 && zqq_m_best < 140") ## loose
    #df = df.Filter("zqq_m_best > 80 && zqq_m_best < 110") ## aggressive 80-60 seems not too much of a difference --> because we fit it!
    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut3"))


    ## CUT ON p(qq)
    hists.append(df.Histo1D(("zqq_p_best_nOne", "", *(200, 0, 200)), "zqq_p_best"))
    ##df = df.Filter("zqq_p_best < 60 && zqq_p_best > 30")
    #df = df.Filter("zqq_p_best < 60") # see significance
    df = df.Filter("zqq_p_best < 90 && zqq_p_best > 20")
    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut4"))


    ## cut on cos(qq)
    hists.append(df.Histo1D(("z_costheta_nOne", "", *bins_cos_abs), "z_costheta"))
    #df = df.Filter("z_costheta < 0.95") # tight cut on pqq
    df = df.Filter("z_costheta < 0.85") # after loose cut on pqq
    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut5"))

    ## acolinearity
    df = df.Define("acolinearity", "FCCAnalyses::acolinearity(zqq_jets_best)")
    hists.append(df.Histo1D(("acolinearity_nOne", "", *bins_aco), "acolinearity"))
    #df = df.Filter("acolinearity > 0.35") # tight cut on pqq
    df = df.Filter("acolinearity > 0.35") # sigificance 0.45, but keep it
    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut6"))


    ## acoplanarity
    df = df.Define("acoplanarity", "FCCAnalyses::acoplanarity(zqq_jets_best)")
    hists.append(df.Histo1D(("acoplanarity_nOne", "", *bins_aco), "acoplanarity"))
    df = df.Filter("acoplanarity < 5")
    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut7"))




    hists.append(df.Histo2D(("W1_m_W2_m_nOne", "", *((150, 0, 150)+(150, 0, 150))), "W1_m", "W2_m"))
    hists.append(df.Histo1D(("delta_mWW_nOne", "", *(1000, 0, 100)), "delta_mWW"))
    hists.append(df.Histo1D(("W1_m_nOne", "", *(200, 0, 200)), "W1_m"))
    hists.append(df.Histo1D(("W2_m_nOne", "", *(200, 0, 200)), "W2_m"))


    df =df.Filter("delta_mWW > 6")
    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut8"))


    hists.append(df.Histo1D(("W1_p_nOne", "", *(200, 0, 200)), "W1_p"))
    hists.append(df.Histo1D(("W2_p_nOne", "", *(200, 0, 200)), "W2_p"))


    ####
    ## CUT 7: cos theta(miss)
    ####
    df = df.Define("missingEnergy_rp", "FCCAnalyses::missingEnergy(ecm, ReconstructedParticles)")
    df = df.Define("missingEnergy", "missingEnergy_rp[0].energy")
    df = df.Define("cosThetaMiss", "FCCAnalyses::get_cosTheta_miss(missingEnergy_rp)")
    hists.append(df.Histo1D(("cosThetaMiss_nOne", "", *bins_cosThetaMiss), "cosThetaMiss"))
    
    df = df.Filter("cosThetaMiss < .995")
    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut9"))

    ## thrust
    df = df.Define("rps_charged_idx", "RP_q != 0") # select only charged particles
    df = df.Define("rps_charged", "ReconstructedParticles[rps_charged_idx]")
    df = df.Define("rps_charged_p", "FCCAnalyses::ReconstructedParticle::get_p(rps_charged)")
    df = df.Define("rps_charged_n", "FCCAnalyses::ReconstructedParticle::get_n(rps_charged)")
    df = df.Define("RP_px_q", "RP_px[rps_charged_idx]")
    df = df.Define("RP_py_q", "RP_py[rps_charged_idx]")
    df = df.Define("RP_pz_q", "RP_pz[rps_charged_idx]")
    df = df.Define("max_p_idx", "int idx=-1; float max=-1; for(int i=0; i<rps_charged_n; i++) if(rps_charged_p[i]>max) {max= rps_charged_p[i]; idx=i;}; return idx")
    
    
    df = df.Define("thrust", "FCCAnalyses::Algorithms::calculate_thrust()(RP_px, RP_py, RP_pz)")
    df = df.Define("thrust_magn", "thrust[0]")
    df = df.Define("thrust_costheta", "abs(thrust[3])")

    
    #df = df.Define("thrust", "FCCAnalyses::minimize_thrust_mc(RP_px_q[max_p_idx], RP_py_q[max_p_idx], RP_pz_q[max_p_idx])(RP_px_q, RP_py_q, RP_pz_q)")
    #df = df.Define("thrust_magn", "thrust[0]")
    #df = df.Define("thrust_costheta", "abs(cos(thrust[5]))")

    hists.append(df.Histo1D(("thrust_magn", "", *(1000, 0, 1)), "thrust_magn"))
    hists.append(df.Histo1D(("thrust_costheta", "", *bins_cos_abs), "thrust_costheta"))

    df = tmva_helper.run_inference(df, col_name="mva_score")
    hists.append(df.Histo1D(("mva_score", "", *(1000, 0, 1)), "mva_score"))

    #df = df.Filter("mva_score[0] > 0.4")
    ## final histograms
    hists.append(df.Histo1D(("zqq_recoil_m", "", *(200, 0, 200)), "zqq_recoil_m_best"))
    hists.append(df.Histo1D(("zqq_m", "", *(200, 0, 200)), "zqq_m_best"))
    #hists.append(df.Histo2D(("zqq_recoil_m_mqq", "", *((50, 100, 150)+(60, 60, 120))), "zqq_recoil_m_best", "zqq_m_best"))
    hists.append(df.Histo3D(("zqq_recoil_m_mqq_pqq", "", *((50, 100, 150)+(100, 40, 140)+(40, 20, 60))), "zqq_recoil_m_best", "zqq_m_best", "zqq_p_best"))
    hists.append(df.Histo2D(("zqq_recoil_m_mqq", "", *((50, 100, 150)+(100, 40, 140))), "zqq_recoil_m_best", "zqq_m_best"))
    hists.append(df.Histo2D(("zqq_recoil_m_pqq", "", *((50, 100, 150)+(40, 20, 60))), "zqq_recoil_m_best", "zqq_p_best"))



if treemaker:
    class RDFanalysis():
        def analysers(df):
            hists, weightsum, df = build_graph_hwwzz(df, "")
            return df

        # define output branches to be saved
        def output():
            branchList = ["zqq_recoil_m_best", "zqq_p_best", "zqq_m_best", "leading_jet_costheta", "subleading_jet_costheta", "leading_jet_p", "subleading_jet_p", "acolinearity", "acoplanarity", "W1_p", "W2_p", "W1_m", "W2_m", "z_costheta", "thrust_magn", "W1_costheta", "W2_costheta"]
            return branchList

else:
    def build_graph(df, dataset):
        hists, weightsum, df = build_graph_hwwzz(df, dataset)
        return hists, weightsum
