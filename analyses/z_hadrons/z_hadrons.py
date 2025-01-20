
import ROOT
ROOT.TH1.SetDefaultSumw2(ROOT.kTRUE)



# list of all processes
fraction = 1
processList = {
    'wz3p6_ee_qq_ecm91p2': {'fraction':fraction},
    'kkmcee_ee_uu_ecm91p2': {'fraction':fraction},
    'kkmcee_ee_dd_ecm91p2': {'fraction':fraction},
    'kkmcee_ee_cc_ecm91p2': {'fraction':fraction},
    'kkmcee_ee_ss_ecm91p2': {'fraction':fraction},
    'kkmcee_ee_bb_ecm91p2': {'fraction':fraction},

    'p8_ee_uu_ecm91p2': {'fraction':fraction},
    'wz3p6_ee_uu_ALEPH_ecm91p2': {'fraction':fraction},
    'wz3p6_ee_uu_ecm91p2': {'fraction':fraction},
    'kkmcee_ee_mumu_ecm91p2': {'fraction':fraction},

    'kkmcee_ee_ee_ecm91p2': {'fraction':fraction},
    'kkmcee_ee_mumu_ecm91p2': {'fraction':fraction},
    'wzp6_ee_mumu_ecm91p2': {'fraction':fraction},

    'p8_ee_gaga_qq_ecm91p2': {'fraction':fraction},
    'wz3p6_ee_gaga_qq_ecm91p2': {'fraction':fraction},

    'wzp6_ee_tautau_ecm91p2': {'fraction':fraction},
    'p8_ee_Ztautau_ecm91': {'fraction':fraction},
    'kkmcee_ee_tautau_ecm91p2': {'fraction':fraction},
}



inputDir = "/ceph/submit/data/group/fcc/ee/generation/DelphesEvents/winter2023/IDEA/"
procDict = "/ceph/submit/data/group/fcc/ee/generation/DelphesEvents/winter2023/IDEA/samplesDict.json"

# additional/custom C++ functions
includePaths = ["../../functions/functions.h", "../../functions/functions_gen.h"]


# output directory
outputDir   = "output/tutorials/z_hadrons/"


# optional: ncpus, default is 4, -1 uses all cores available
nCPUS       = 256

# scale the histograms with the cross-section and integrated luminosity
doScale = True
intLumi = 100e6 # 100 ab-1

# define histograms
bins_p_mu = (2000, 0, 200) # 100 MeV bins
bins_m_ll = (2000, 0, 200) # 100 MeV bins
bins_p_ll = (2000, 0, 200) # 100 MeV bins

bins_theta = (315, 0, 3.15)
bins_cos_abs = (100, 0, 1)
bins_phi = (500, -5, 5)

bins_count = (100, 0, 100)
bins_pdgid = (60, -30, 30)
bins_charge = (10, -5, 5)

dijet_m = (200, 0, 200) # 1 GeV bins

bins_p = (1000, 0, 100)
bins_m = (100, 0, 100)
bins_norm = (200, 0, 2)
bins_nparticles = (200, 0, 200)



def build_graph(df, dataset):

    hists = []
    df = df.Define("ecm", "91.2")

    df = df.Define("weight", "1.0")
    weightsum = df.Sum("weight")


    df = df.Alias("Particle0", "Particle#0.index")
    df = df.Alias("Particle1", "Particle#1.index")
    df = df.Alias("MCRecoAssociations0", "MCRecoAssociations#0.index")
    df = df.Alias("MCRecoAssociations1", "MCRecoAssociations#1.index")

    # get reco particles
    df = df.Define("RP_px", "FCCAnalyses::ReconstructedParticle::get_px(ReconstructedParticles)")
    df = df.Define("RP_py", "FCCAnalyses::ReconstructedParticle::get_py(ReconstructedParticles)")
    df = df.Define("RP_pz", "FCCAnalyses::ReconstructedParticle::get_pz(ReconstructedParticles)")
    df = df.Define("RP_p", "FCCAnalyses::ReconstructedParticle::get_p(ReconstructedParticles)")
    df = df.Define("RP_e",  "FCCAnalyses::ReconstructedParticle::get_e(ReconstructedParticles)")
    df = df.Define("RP_m",  "FCCAnalyses::ReconstructedParticle::get_mass(ReconstructedParticles)")
    df = df.Define("RP_q",  "FCCAnalyses::ReconstructedParticle::get_charge(ReconstructedParticles)")
    df = df.Define("RP_theta", "FCCAnalyses::ReconstructedParticle::get_theta(ReconstructedParticles)")
    df = df.Define("RP_no",  "FCCAnalyses::ReconstructedParticle::get_n(ReconstructedParticles)")

    df = df.Define("ReconstructedParticles_q", "FCCAnalyses::ReconstructedParticle::sel_charge(1, true)(ReconstructedParticles)")
    df = df.Define("RP_q_p", "FCCAnalyses::ReconstructedParticle::get_p(ReconstructedParticles_q)")
    df = df.Define("RP_q_no",  "FCCAnalyses::ReconstructedParticle::get_n(ReconstructedParticles_q)")

    # get gen particles
    df = df.Define("ParticleStable", "FCCAnalyses::MCParticle::sel_genStatus(1)(Particle)")
    df = df.Define("GP_px", "FCCAnalyses::MCParticle::get_px(ParticleStable)")
    df = df.Define("GP_py", "FCCAnalyses::MCParticle::get_py(ParticleStable)")
    df = df.Define("GP_pz", "FCCAnalyses::MCParticle::get_pz(ParticleStable)")
    df = df.Define("GP_p", "FCCAnalyses::MCParticle::get_p(ParticleStable)")
    df = df.Define("GP_e",  "FCCAnalyses::MCParticle::get_e(ParticleStable)")
    df = df.Define("GP_m",  "FCCAnalyses::MCParticle::get_mass(ParticleStable)")
    df = df.Define("GP_q",  "FCCAnalyses::MCParticle::get_charge(ParticleStable)")
    df = df.Define("GP_theta",  "FCCAnalyses::MCParticle::get_theta(ParticleStable)")
    df = df.Define("GP_no",  "FCCAnalyses::MCParticle::get_n(ParticleStable)")
    df = df.Define("Particle_q", "FCCAnalyses::Vec_mc ret; for(int i=0; i<GP_no; i++) if(GP_q[i]!=0) ret.push_back(ParticleStable[i]); return ret;")
    df = df.Define("GP_q_p", "FCCAnalyses::MCParticle::get_p(Particle_q)")
    df = df.Define("GP_q_no",  "FCCAnalyses::MCParticle::get_n(Particle_q)")
    
    df = df.Define("visibleMass", "FCCAnalyses::visibleMass(ReconstructedParticles)")
    df = df.Define("missingEnergy_vec", "FCCAnalyses::missingEnergy(ecm, ReconstructedParticles)")
    df = df.Define("missingEnergy", "FCCAnalyses::ReconstructedParticle::get_e(missingEnergy_vec)")
    df = df.Define("visible_energy", "FCCAnalyses::visibleEnergy(ReconstructedParticles)")
    df = df.Define("visible_energy_norm", "visible_energy/ecm")


    # calculate thrust axis (based on charged particles)
    df = df.Define("max_p_idx", "int idx=-1; float max=-1; for(int i=0; i<RP_q_no; i++) if(RP_q_p[i]>max) {max=RP_q_p[i]; idx=i;}; return idx")

    #df = df.Define("thrust", "FCCAnalyses::minimize_thrust_mc()(RP_px, RP_py, RP_pz)")
    #df = df.Define("thrust", "FCCAnalyses::minimize_thrust_mc(RP_px_q[max_p_idx], RP_py_q[max_p_idx], RP_pz_q[max_p_idx])(RP_px_q, RP_py_q, RP_pz_q)")
    #df = df.Define("thrust", "FCCAnalyses::minimize_thrust_mc()(RP_px_q, RP_py_q, RP_pz_q)")
    #df = df.Define("thrust", "FCCAnalyses::minimize_thrust_mc(1, 1, 1)(RP_px, RP_py, RP_pz)")
    #df = df.Define("thrust_magn", "thrust[0]")
    #df = df.Define("thrust_costheta", "abs(cos(thrust[5]))")

    # second way to calculate thrust
    df = df.Define("thrust", "FCCAnalyses::Algorithms::calculate_thrust()(RP_px, RP_py, RP_pz)")
    df = df.Define("thrust_magn", "thrust[0]")
    df = df.Define("thrust_costheta", "abs(thrust[3])")

    hists.append(df.Histo1D(("thrust_magn_nCut", "", *bins_norm), "thrust_magn"))
    hists.append(df.Histo1D(("thrust_costheta_nCut", "", *bins_cos_abs), "thrust_costheta"))

    df = df.Define("energy_imbalance", "FCCAnalyses::energy_imbalance(ReconstructedParticles)")
    df = df.Define("energy_imbalance_tot", "energy_imbalance[0]")
    df = df.Define("energy_imbalance_trans", "energy_imbalance[1]/energy_imbalance[0]")
    df = df.Define("energy_imbalance_long", "energy_imbalance[2]/energy_imbalance[0]")


    # jet clustering
    # more info: https://indico.cern.ch/event/1173562/contributions/4929025/attachments/2470068/4237859/2022-06-FCC-jets.pdf
    # https://github.com/HEP-FCC/FCCAnalyses/blob/master/addons/FastJet/src/JetClustering.cc
    df = df.Define("pseudo_jets", "FCCAnalyses::JetClusteringUtils::set_pseudoJets(RP_px, RP_py, RP_pz, RP_e)")
    df = df.Define("clustered_jets", "JetClustering::clustering_ee_kt(2, 2, 1, 0)(pseudo_jets)")
    df = df.Define("jets", "FCCAnalyses::JetClusteringUtils::get_pseudoJets(clustered_jets)")
    df = df.Define("jetconstituents", "FCCAnalyses::JetClusteringUtils::get_constituents(clustered_jets)") # one-to-one mapping to reconstructedparticles
    df = df.Define("jets_e", "FCCAnalyses::JetClusteringUtils::get_e(jets)")
    df = df.Define("jets_px", "FCCAnalyses::JetClusteringUtils::get_px(jets)")
    df = df.Define("jets_py", "FCCAnalyses::JetClusteringUtils::get_py(jets)")
    df = df.Define("jets_pz", "FCCAnalyses::JetClusteringUtils::get_pz(jets)")
    df = df.Define("jets_phi", "FCCAnalyses::JetClusteringUtils::get_phi(jets)")
    df = df.Define("jets_m", "FCCAnalyses::JetClusteringUtils::get_m(jets)")

    df = df.Define("jets_tlv", "FCCAnalyses::makeLorentzVectors(jets_px, jets_py, jets_pz, jets_e)")
    df = df.Define("njets", "jets_e.size()")
    df = df.Define("zqq", "jets_tlv[0]+jets_tlv[1]")
    df = df.Define("zqq_m", "zqq.M()")

    hists.append(df.Histo1D(("njets", "", *bins_count), "njets"))
    hists.append(df.Histo1D(("zqq_m", "", *bins_m), "zqq_m"))

    hists.append(df.Histo1D(("RP_no", "", *bins_count), "RP_no"))
    hists.append(df.Histo1D(("GP_no", "", *bins_count), "GP_no"))
    hists.append(df.Histo1D(("RP_q_no", "", *bins_count), "RP_q_no"))
    hists.append(df.Histo1D(("GP_q_no", "", *bins_count), "GP_q_no"))

    hists.append(df.Histo1D(("RP_p", "", *bins_p), "RP_p"))
    hists.append(df.Histo1D(("GP_p", "", *bins_p), "GP_p"))
    hists.append(df.Histo1D(("RP_q_p", "", *bins_p), "RP_q_p"))
    hists.append(df.Histo1D(("GP_q_p", "", *bins_p), "GP_q_p"))

    hists.append(df.Histo1D(("RP_theta", "", *bins_theta), "RP_theta"))
    hists.append(df.Histo1D(("GP_theta", "", *bins_theta), "GP_theta"))

    df = df.Define("cut0", "0")
    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut0"))


    ####################
    # plots without cuts
    df_b = df.Filter("thrust_costheta < 0.74")
    df_e = df.Filter("thrust_costheta > 0.74")
    hists.append(df_b.Histo1D(("RP_no_barrel_nCut", "", *bins_nparticles), "RP_no"))
    hists.append(df_e.Histo1D(("RP_no_endcap_nCut", "", *bins_nparticles), "RP_no"))
    hists.append(df.Histo1D(("visible_energy_norm_nCut", "", *bins_norm), "visible_energy_norm"))
    hists.append(df.Histo1D(("energy_imbalance_long_nCut", "", *bins_norm), "energy_imbalance_long"))
    hists.append(df.Histo1D(("energy_imbalance_trans_nCut", "", *bins_norm), "energy_imbalance_trans"))


    ####################
    # N-1 plots
    sel_visible_energy = "(visible_energy_norm < 2 && visible_energy_norm > 0.5)"
    sel_energy_imbalance_long = "(energy_imbalance_long < 0.6)"
    sel_energy_imbalance_trans = "(energy_imbalance_trans < 0.6)"
    sel_nparticles_barrel = "(RP_no > 13)"
    sel_nparticles_endcap = "(RP_no > 17)"
    sel_barrel = "(thrust_costheta < 0.74)"
    sel_endcap = "(thrust_costheta > 0.74)"
    sel_nparticles = f"(({sel_barrel} && {sel_nparticles_barrel}) || ({sel_endcap} && {sel_nparticles_endcap}))"

    df_nparticles_barrel = df.Filter(f"{sel_barrel} && {sel_visible_energy} && {sel_energy_imbalance_long} && {sel_energy_imbalance_trans}")
    df_nparticles_endcap = df.Filter(f"{sel_endcap} && {sel_visible_energy} && {sel_energy_imbalance_long} && {sel_energy_imbalance_trans}")
    hists.append(df_nparticles_barrel.Histo1D(("RP_no_barrel_nOne", "", *bins_nparticles), "RP_no"))
    hists.append(df_nparticles_endcap.Histo1D(("RP_no_endcap_nOne", "", *bins_nparticles), "RP_no"))

    df_visible_energy = df.Filter(f"{sel_energy_imbalance_long} && {sel_energy_imbalance_trans} && {sel_nparticles}")
    hists.append(df_visible_energy.Histo1D(("visible_energy_norm_nOne", "", *bins_norm), "visible_energy_norm"))

    df_energy_imbalance_long = df.Filter(f"{sel_visible_energy} && {sel_energy_imbalance_trans} && {sel_nparticles}")
    hists.append(df_energy_imbalance_long.Histo1D(("energy_imbalance_long_nOne", "", *bins_norm), "energy_imbalance_long"))

    df_energy_imbalance_trans = df.Filter(f"{sel_visible_energy} && {sel_energy_imbalance_long} && {sel_nparticles}")
    hists.append(df_energy_imbalance_trans.Histo1D(("energy_imbalance_trans_nOne", "", *bins_norm), "energy_imbalance_trans"))



    # all filters
    df = df.Filter(f"{sel_visible_energy} && {sel_energy_imbalance_long} && {sel_energy_imbalance_trans} && {sel_nparticles}")

    hists.append(df.Histo1D(("visible_energy_norm", "", *bins_norm), "visible_energy_norm"))
    hists.append(df.Histo1D(("energy_imbalance_trans", "", *bins_norm), "energy_imbalance_trans"))
    hists.append(df.Histo1D(("energy_imbalance_long", "", *bins_norm), "energy_imbalance_long"))
    
    hists.append(df.Histo1D(("thrust_magn", "", *bins_norm), "thrust_magn"))
    hists.append(df.Histo1D(("thrust_costheta", "", *bins_cos_abs), "thrust_costheta"))

    hists.append(df.Histo1D(("visibleMass", "", *bins_m), "visibleMass"))
    hists.append(df.Histo1D(("missingEnergy", "", *bins_m), "missingEnergy"))

    df = df.Define("cut1", "1")
    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut1"))

    return hists, weightsum
