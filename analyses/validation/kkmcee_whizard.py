
import ROOT
ROOT.TH1.SetDefaultSumw2(ROOT.kTRUE)



# list of all processes
fraction = 0.05
processList = {
    'kkmcee_ee_uu_ecm240':          {'fraction':fraction},
    'kkmcee_ee_mumu_ecm240':        {'fraction':fraction},
    'wz3p6_ee_uu_ecm240':           {'fraction':fraction},
    'wzp6_ee_mumu_ecm240':          {'fraction':fraction},
}





inputDir = "/ceph/submit/data/group/fcc/ee/generation/DelphesEvents/winter2023/IDEA/"
procDict = "/ceph/submit/data/group/fcc/ee/generation/DelphesEvents/winter2023/IDEA/samplesDict.json"

# additional/custom C++ functions
includePaths = ["../../functions/functions.h", "../../functions/functions_gen.h"]


# output directory
outputDir   = "output/kkmcee_whizard/histmaker/"


# optional: ncpus, default is 4, -1 uses all cores available
nCPUS       = 32

# scale the histograms with the cross-section and integrated luminosity
doScale = True
intLumi = 10800000

# define histograms
bins_m = (250, 0, 250) # 100 MeV bins
bins_maa = (100, 120, 130) # 100 MeV bins
bins_p_photons = (20000, 0, 20) # 100 MeV bins
bins_p = (2500, 0, 250) # 100 MeV bins
bins_m_zoom = (100, 120, 130) # 100 MeV
bins_m_ll = (25000, 0, 250) # 10 MeV bins
bins_p_ll = (25000, 0, 250) # 10 MeV bins

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

bins_reso = (10000, 0.95, 1.05)
bins_ecm_eff = (1000, 0.0, 1.0)
bins_cos = (100, -1, 1)
bins_acolinearity_deg = (1000, 0.0, 90.0)
bins_acolinearity_rad = (1000, 0.0, 1.0)



def build_graph(df, dataset):

    hists, cols = [], []

    df = df.Define("weight", "1.0")
    weightsum = df.Sum("weight")



    # define collections
    df = df.Alias("Particle0", "Particle#0.index")
    df = df.Alias("Particle1", "Particle#1.index")
    df = df.Alias("MCRecoAssociations0", "MCRecoAssociations#0.index")
    df = df.Alias("MCRecoAssociations1", "MCRecoAssociations#1.index")
    df = df.Alias("Muon0", "Muon#0.index")
    df = df.Alias("Photon0", "Photon#0.index")

    # gen photons
    df = df.Define("gen_photons", "FCCAnalyses::get_gen_pdg(Particle, 22)")
    df = df.Define("gen_photons_p", "FCCAnalyses::MCParticle::get_p(gen_photons)")
    df = df.Define("gen_photons_theta", "FCCAnalyses::MCParticle::get_theta(gen_photons)")
    df = df.Define("gen_photons_phi", "FCCAnalyses::MCParticle::get_phi(gen_photons)")
    df = df.Define("gen_photons_no", "FCCAnalyses::MCParticle::get_n(gen_photons)")

    hists.append(df.Histo1D(("gen_photons_p", "", *bins_p_photons), "gen_photons_p"))
    hists.append(df.Histo1D(("gen_photons_theta", "", *bins_theta), "gen_photons_theta"))
    hists.append(df.Histo1D(("gen_photons_phi", "", *bins_phi), "gen_photons_phi"))
    hists.append(df.Histo1D(("gen_photons_no", "", *bins_count), "gen_photons_no"))

    df = df.Define("gen_photon1_p", "gen_photons_p[0]")
    df = df.Define("gen_photon2_p", "gen_photons_p[1]")
    
    hists.append(df.Histo1D(("gen_photon1_p", "", *bins_p), "gen_photon1_p"))
    hists.append(df.Histo1D(("gen_photon2_p", "", *bins_p), "gen_photon2_p"))


    # reco photons
    df = df.Define("photons", "FCCAnalyses::ReconstructedParticle::get(Photon0, ReconstructedParticles)")
    df = df.Define("photons_p", "FCCAnalyses::ReconstructedParticle::get_p(photons)")
    df = df.Define("photons_theta", "FCCAnalyses::ReconstructedParticle::get_theta(photons)")
    df = df.Define("photons_phi", "FCCAnalyses::ReconstructedParticle::get_phi(photons)")
    df = df.Define("photons_no", "FCCAnalyses::ReconstructedParticle::get_n(photons)")
    hists.append(df.Histo1D(("photons_p", "", *bins_p), "photons_p"))
    hists.append(df.Histo1D(("photons_theta", "", *bins_theta), "photons_theta"))
    hists.append(df.Histo1D(("photons_phi", "", *bins_phi), "photons_phi"))
    hists.append(df.Histo1D(("photons_no", "", *bins_count), "photons_no"))

    df = df.Define("photon1_p", "photons_p[0]")
    df = df.Define("photon2_p", "photons_p[1]")
    hists.append(df.Histo1D(("photon1_p", "", *bins_p), "photon1_p"))
    hists.append(df.Histo1D(("photon2_p", "", *bins_p), "photon2_p"))

    # reco leptons
    df = df.Define("leptons", "FCCAnalyses::ReconstructedParticle::get(Muon0, ReconstructedParticles)")
    df = df.Define("leptons_p", "FCCAnalyses::ReconstructedParticle::get_p(leptons)")
    df = df.Define("leptons_theta", "FCCAnalyses::ReconstructedParticle::get_theta(leptons)")
    df = df.Define("leptons_phi", "FCCAnalyses::ReconstructedParticle::get_phi(leptons)")
    df = df.Define("leptons_charge", "FCCAnalyses::ReconstructedParticle::get_charge(leptons)")
    df = df.Define("leptons_no", "FCCAnalyses::ReconstructedParticle::get_n(leptons)")
    hists.append(df.Histo1D(("leptons_p", "", *bins_p), "leptons_p"))
    hists.append(df.Histo1D(("leptons_theta", "", *bins_theta), "leptons_theta"))
    hists.append(df.Histo1D(("leptons_phi", "", *bins_phi), "leptons_phi"))
    hists.append(df.Histo1D(("leptons_charge", "", *bins_charge), "leptons_charge"))
    hists.append(df.Histo1D(("leptons_no", "", *bins_count), "leptons_no"))


    
    # reconstruct resonance
    df = df.Filter("leptons_no == 2")
    df = df.Define("leptons_tlv", "FCCAnalyses::makeLorentzVectors(leptons)")
    df = df.Define("m_ll", "(leptons_tlv[0]+leptons_tlv[1]).M()")
    df = df.Define("p_ll", "(leptons_tlv[0]+leptons_tlv[1]).P()")
    hists.append(df.Histo1D(("m_ll", "", *bins_m_ll), "m_ll"))
    hists.append(df.Histo1D(("p_ll", "", *bins_p_ll), "p_ll"))
  

    # muon angular analysis (effective ECM, AFB, acolinearity, ...)
    df = df.Define("theta_plus", "(leptons_charge[0] > 0) ? leptons_theta[0] : leptons_theta[1]")
    df = df.Define("theta_minus", "(leptons_charge[0] < 0) ? leptons_theta[0] : leptons_theta[1]")
    df = df.Define("cos_theta_plus", "cos(theta_plus)")
    df = df.Define("cos_theta_minus", "cos(theta_minus)")
    df = df.Define("acolinearity_rad", "FCCAnalyses::acolinearity(leptons)")
    df = df.Define("acolinearity_deg", "acolinearity_rad*180./acos(-1)")
    
    df = df.Define("cosThetac", "(sin(theta_plus-theta_minus))/(sin(theta_plus)+sin(theta_minus))")
    df = df.Define("ecm_eff", "(sin(theta_plus) + sin(theta_minus) - sin(theta_plus+theta_minus))/(sin(theta_plus) + sin(theta_minus) + sin(theta_plus+theta_minus))")
    
    hists.append(df.Histo1D(("theta_plus", "", *bins_theta), "theta_plus"))
    hists.append(df.Histo1D(("theta_minus", "", *bins_theta), "theta_minus"))
    hists.append(df.Histo1D(("cos_theta_plus", "", *bins_cos), "cos_theta_plus"))
    hists.append(df.Histo1D(("cos_theta_minus", "", *bins_cos), "cos_theta_minus"))
    hists.append(df.Histo1D(("cosThetac", "", *bins_cos), "cosThetac"))
    hists.append(df.Histo1D(("ecm_eff", "", *bins_ecm_eff), "ecm_eff"))
    hists.append(df.Histo1D(("acolinearity_rad", "", *bins_acolinearity_rad), "acolinearity_rad"))
    hists.append(df.Histo1D(("acolinearity_deg", "", *bins_acolinearity_deg), "acolinearity_deg"))



    
    
    df = df.Define("rps_no_muons", "FCCAnalyses::ReconstructedParticle::remove(ReconstructedParticles, leptons)")
    df = df.Define("visibleEnergy", "FCCAnalyses::visibleEnergy(rps_no_muons)")
    hists.append(df.Histo1D(("visibleEnergy", "", *bins_p_ll), "visibleEnergy"))

    # radiative return events
    df = df.Filter("m_ll < 93 && m_ll > 90")
    #df = df.Filter("m_ll > 220")
    
    hists.append(df.Histo1D(("p_ll_rr", "", *bins_p_ll), "p_ll"))
    hists.append(df.Histo1D(("visibleEnergy_rr", "", *bins_p_ll), "visibleEnergy"))
    
    df = df.Define("leading", "(photon1_p > photon2_p) ? photon1_p : photon2_p")
    df = df.Define("subleading", "(photon1_p > photon2_p) ? photon2_p : photon1_p")
    hists.append(df.Histo1D(("photon1_p_rr", "", *bins_p), "leading"))
    hists.append(df.Histo1D(("photon2_p_rr", "", *bins_p), "subleading"))
    


    return hists, weightsum
