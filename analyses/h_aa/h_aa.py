
import ROOT
ROOT.TH1.SetDefaultSumw2(ROOT.kTRUE)

ecm = 240

# list of all processes
fraction = 0.01
fraction = 1


if ecm == 240:
    processList = {
        'wzp6_ee_nunuH_Haa_ecm240':         {'fraction':1},
        'wzp6_ee_eeH_Haa_ecm240':           {'fraction':1},
        'wzp6_ee_tautauH_Haa_ecm240':       {'fraction':1},
        'wzp6_ee_ccH_Haa_ecm240':           {'fraction':1},
        'wzp6_ee_bbH_Haa_ecm240':           {'fraction':1},
        'wzp6_ee_qqH_Haa_ecm240':           {'fraction':1},
        'wzp6_ee_ssH_Haa_ecm240':           {'fraction':1},
        'wzp6_ee_mumuH_Haa_ecm240':         {'fraction':1},
        'wzp6_ee_nuenueH_Haa_ecm240':       {'fraction':fraction},
        'wzp6_ee_numunumuH_Haa_ecm240':     {'fraction':fraction},
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
    }
    processList = {
        'wzp6_ee_nunuH_Haa_ecm240':         {'fraction':1},
        'wzp6_ee_eeH_Haa_ecm240':           {'fraction':1},
        'wzp6_ee_tautauH_Haa_ecm240':       {'fraction':1},
        'wzp6_ee_ccH_Haa_ecm240':           {'fraction':1},
        'wzp6_ee_bbH_Haa_ecm240':           {'fraction':1},
        'wzp6_ee_qqH_Haa_ecm240':           {'fraction':1},
        'wzp6_ee_ssH_Haa_ecm240':           {'fraction':1},
        'wzp6_ee_mumuH_Haa_ecm240':         {'fraction':1},
        'wzp6_ee_nuenueH_Haa_ecm240':       {'fraction':fraction},
        'wzp6_ee_numunumuH_Haa_ecm240':     {'fraction':fraction},
    }

if ecm == 365:
    processList = {
        'wzp6_ee_nunuH_Haa_ecm365':         {'fraction':1},
        'wzp6_ee_eeH_Haa_ecm365':           {'fraction':1},
        'wzp6_ee_tautauH_Haa_ecm365':       {'fraction':1},
        'wzp6_ee_ccH_Haa_ecm365':           {'fraction':1},
        'wzp6_ee_bbH_Haa_ecm365':           {'fraction':1},
        'wzp6_ee_qqH_Haa_ecm365':           {'fraction':1},
        'wzp6_ee_ssH_Haa_ecm365':           {'fraction':1},
        'wzp6_ee_mumuH_Haa_ecm365':         {'fraction':1},
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
        'wzp6_ee_nuenueH_Haa_ecm365':       {'fraction':fraction},
        'wzp6_ee_numunumuH_Haa_ecm365':     {'fraction':fraction},
    }

    processList = {
        'wzp6_ee_nunuH_Haa_ecm365':         {'fraction':1},
        'wzp6_ee_eeH_Haa_ecm365':           {'fraction':1},
        'wzp6_ee_tautauH_Haa_ecm365':       {'fraction':1},
        'wzp6_ee_ccH_Haa_ecm365':           {'fraction':1},
        'wzp6_ee_bbH_Haa_ecm365':           {'fraction':1},
        'wzp6_ee_qqH_Haa_ecm365':           {'fraction':1},
        'wzp6_ee_ssH_Haa_ecm365':           {'fraction':1},
        'wzp6_ee_mumuH_Haa_ecm365':         {'fraction':1},
        'wzp6_ee_nuenueH_Haa_ecm365':       {'fraction':fraction},
        'wzp6_ee_numunumuH_Haa_ecm365':     {'fraction':fraction},
    }




inputDir = "/ceph/submit/data/group/fcc/ee/generation/DelphesEvents/winter2023/IDEA/"
procDict = "/ceph/submit/data/group/fcc/ee/generation/DelphesEvents/winter2023/IDEA/samplesDict.json"

# additional/custom C++ functions
includePaths = ["../../functions/functions.h", "../../functions/functions_gen.h"]


# output directory
outputDir   = f"output/h_aa/histmaker/ecm{ecm}/"

nCPUS       = 8

# scale the histograms with the cross-section and integrated luminosity
doScale = True
intLumi = 10800000 if ecm == 240 else 3e6

# define histograms
bins_m = (500, 0, 500) # 100 MeV bins
bins_maa = (100, 120, 130) # 100 MeV bins
bins_p = (500, 0, 500) # 100 MeV bins
bins_m_zoom = (100, 120, 130) # 100 MeV

bins_theta = (500, 0, 5)
bins_phi = (400, -4, 4)

bins_count = (100, 0, 100)
bins_pdgid = (60, -30, 30)
bins_charge = (10, -5, 5)

bins_resolution = (10000, 0.95, 1.05)
bins_resolution_1 = (20000, 0, 2)

jet_energy = (5000, 0, 500) # 100 MeV bins
dijet_m = (5000, 0, 500) # 100 MeV bins
visMass = (5000, 0, 500) # 100 MeV bins
missEnergy  = (5000, 0, 500) # 100 MeV bins

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

    df = df.Define("ecm", f"{ecm}")

    '''
    Some events have ReconstructedParticles.size  > MCRecoAssociations.size. In these events there is an extra RECO photon that might be added by Delphes. It's a very soft photon that will not affect the physics studies and can be safely neglected. I checked in a few samples and the probability such events happen is < 0.1%. So either we filter out these events or you find a hack in the MC-RECO particle association on how to deal with it 
    '''
    #df = df.Filter("ReconstructedParticles.size() == MCRecoAssociations0.size()")
    #df = df.Filter("if(ReconstructedParticles.size() != MCRecoAssociations0.size()) {std::cout << \"NOT_EQUAL\" << std::endl;}; return (ReconstructedParticles.size() == MCRecoAssociations0.size())")
    #df = df.Filter("MCRecoAssociations0.size() == MCRecoAssociations1.size()")


    # all photons
    df = df.Alias("Photon0", "Photon#0.index")
    df = df.Define("photons_all", "FCCAnalyses::ReconstructedParticle::get(Photon0, ReconstructedParticles)")
    df = df.Define("photons_all_p", "FCCAnalyses::ReconstructedParticle::get_p(photons_all)")
    df = df.Define("photons_all_theta", "FCCAnalyses::ReconstructedParticle::get_theta(photons_all)")
    df = df.Define("photons_all_no", "FCCAnalyses::ReconstructedParticle::get_n(photons_all)")

    hists.append(df.Histo1D(("photons_all_p", "", *bins_p), "photons_all_p"))



    # select photon momentum 40-95 GeV
    if ecm == 240:
        df = df.Define("photons", "FCCAnalyses::sel_range(40, 95, false)(photons_all, photons_all_p)")
    if ecm == 365:
        df = df.Define("photons", "FCCAnalyses::sel_range(20, 170, false)(photons_all, photons_all_p)")
    df = df.Define("photons_p", "FCCAnalyses::ReconstructedParticle::get_p(photons)")
    df = df.Define("photons_theta", "FCCAnalyses::ReconstructedParticle::get_theta(photons)")
    df = df.Define("photons_n", "FCCAnalyses::ReconstructedParticle::get_n(photons)")
    df = df.Define("photons_orig", "photons")

    # do photon smearing
    if "_Haa_" in dataset:
        df = df.Define("photons_idx", "FCCAnalyses::getIndex(photons, ReconstructedParticles)")
        #df = df.Define("photons_smeared_p", "FCCAnalyses::smearPhotonEnergyResolution(photons, photons_idx, MCRecoAssociations0, MCRecoAssociations1, ReconstructedParticles, Particle)")
        #df = df.Redefine("photons", "photons_smeared_p")

        ## compute the photon energy resolution
        df = df.Define("photons_reso_p", "FCCAnalyses::photonResolution(photons, photons_idx, MCRecoAssociations0, MCRecoAssociations1, ReconstructedParticles, Particle, 0)")
        df = df.Define("photons_reso_theta", "FCCAnalyses::photonResolution(photons, photons_idx, MCRecoAssociations0, MCRecoAssociations1, ReconstructedParticles, Particle, 1)")
        df = df.Define("photons_reso_phi", "FCCAnalyses::photonResolution(photons, photons_idx, MCRecoAssociations0, MCRecoAssociations1, ReconstructedParticles, Particle, 2)")
        hists.append(df.Histo1D(("photons_reso_p", "", *bins_resolution), "photons_reso_p"))
        hists.append(df.Histo1D(("photons_reso_theta", "", *bins_resolution), "photons_reso_theta"))
        hists.append(df.Histo1D(("photons_reso_phi", "", *bins_resolution), "photons_reso_phi"))





    #########
    ### CUT 0: all events
    #########
    hists.append(df.Histo1D(("cutFlow_zh", "", *bins_count), "cut0"))
    hists.append(df.Histo1D(("cutFlow_vbf", "", *bins_count), "cut0"))
    hists.append(df.Histo1D(("photons_n", "", *bins_count), "photons_n"))
    hists.append(df.Histo1D(("photons_p", "", *bins_p), "photons_p"))
    

    #########
    ### CUT 1: at least 2 photons, and select the Higgs candidate photons
    #########
    hists.append(df.Histo1D(("photons_n_nOne", "", *bins_count), "photons_n"))
    df = df.Filter("photons_n >= 2")

    # build the Z resonance based on the available leptons. Returns the best lepton pair compatible with the Z mass and recoil at 125 GeV
    # technically, it returns a ReconstructedParticleData object with index 0 the di-lepton system, index and 2 the leptons of the pair
    df = df.Define("hbuilder_result", "FCCAnalyses::resonanceBuilder_mass_recoil(125, 91.2, 0.5, ecm, false)(photons, MCRecoAssociations0, MCRecoAssociations1, ReconstructedParticles, Particle, Particle0, Particle1)")
    df = df.Define("haa_cand", "ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>{hbuilder_result[0]}") # the H
    df = df.Define("haa_photons", "ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>{hbuilder_result[1],hbuilder_result[2]}") # the photons
    df = df.Define("haa_m", "FCCAnalyses::ReconstructedParticle::get_mass(haa_cand)[0]")
    df = df.Define("haa_p", "FCCAnalyses::ReconstructedParticle::get_p(haa_cand)[0]")
    df = df.Define("haa_theta", "FCCAnalyses::ReconstructedParticle::get_theta(haa_cand)[0]")
    df = df.Define("haa_recoil", "FCCAnalyses::ReconstructedParticle::recoilBuilder(ecm)(haa_photons)")
    df = df.Define("haa_recoil_m", "FCCAnalyses::ReconstructedParticle::get_mass(haa_recoil)[0]")


    df = df.Define("haa_photons_p", "FCCAnalyses::ReconstructedParticle::get_p(haa_photons)")
    df = df.Define("haa_photons_theta", "FCCAnalyses::ReconstructedParticle::get_theta(haa_photons)")
    df = df.Define("haa_photons_costheta", "Vec_f ret; for(auto & theta: haa_photons_theta) ret.push_back(std::abs(cos(theta))); return ret;")
    df = df.Define("haa_photons_iso", "FCCAnalyses::coneIsolation(0.01, 0.5)(haa_photons, ReconstructedParticles)")

    df = df.Define("photon_leading_idx", "(haa_photons_p[0] > haa_photons_p[1]) ? 0 : 1")
    df = df.Define("photon_leading_p", "haa_photons_p[photon_leading_idx]")
    df = df.Define("photon_leading_costheta", "haa_photons_costheta[photon_leading_idx]")
    df = df.Define("photon_leading_iso", "haa_photons_iso[photon_leading_idx]")
    
    df = df.Define("photon_subleading_idx", "(haa_photons_p[0] > haa_photons_p[1]) ? 1 : 0")
    df = df.Define("photon_subleading_p", "haa_photons_p[photon_subleading_idx]")
    df = df.Define("photon_subleading_costheta", "haa_photons_costheta[photon_subleading_idx]")
    df = df.Define("photon_subleading_iso", "haa_photons_iso[photon_subleading_idx]")

    hists.append(df.Histo1D(("photon_leading_p", "", *bins_p), "photon_leading_p"))
    hists.append(df.Histo1D(("photon_leading_costheta", "", *bins_cos_abs), "photon_leading_costheta"))
    hists.append(df.Histo1D(("photon_leading_iso", "", *bins_iso), "photon_leading_iso"))


    # do some smearing

    #########
    ### CUT 2: properties of the photons (photon ID)
    #########
    hists.append(df.Histo1D(("cutFlow_zh", "", *bins_count), "cut1"))
    hists.append(df.Histo1D(("cutFlow_vbf", "", *bins_count), "cut1"))
    if ecm == 240:
        #df = df.Filter("photon_leading_iso < 0.1 && photon_leading_costheta < 0.75")
        df = df.Filter("photon_leading_iso < 0.3 && photon_leading_costheta < 0.85")
    if ecm == 365:
        df = df.Filter("photon_leading_iso < 0.1 && photon_leading_costheta < 0.85")


    hists.append(df.Histo1D(("photon_subleading_p", "", *bins_p), "photon_subleading_p"))
    hists.append(df.Histo1D(("photon_subleading_costheta", "", *bins_cos_abs), "photon_subleading_costheta"))
    hists.append(df.Histo1D(("photon_subleading_iso", "", *bins_iso), "photon_subleading_iso"))

    if ecm == 240:
        df = df.Filter("photon_subleading_costheta < 0.85")
    if ecm == 365:
        df = df.Filter("photon_subleading_costheta < 0.95")
    hists.append(df.Histo1D(("cutFlow_zh", "", *bins_count), "cut2"))
    hists.append(df.Histo1D(("cutFlow_vbf", "", *bins_count), "cut2"))


    # missing mass to split ZH/VBF at 365 GeV
    # no distinction at 240 GeV
    df = df.Define("missingMass", "FCCAnalyses::missingMass(ecm, ReconstructedParticles)")
    hists.append(df.Histo1D(("missingMass", "", *bins_m), "missingMass"))
    hists.append(df.Histo1D(("haa_recoil_m_nOne", "", *bins_m), "haa_recoil_m"))
    
    ## thrust
    #df = df.Define("RP_px_all", "FCCAnalyses::ReconstructedParticle::get_px(ReconstructedParticles)")
    #df = df.Define("RP_py_all", "FCCAnalyses::ReconstructedParticle::get_py(ReconstructedParticles)")
    #df = df.Define("RP_pz_all","FCCAnalyses::ReconstructedParticle::get_pz(ReconstructedParticles)")
    #df = df.Define("thrust", "FCCAnalyses::Algorithms::calculate_thrust()(RP_px_all, RP_py_all, RP_pz_all)")
    #df = df.Define("thrust_magn", "thrust[0]")
    #df = df.Define("thrust_costheta", "abs(thrust[3])")


    #hists.append(df.Histo1D(("thrust_magn", "", *(1000, 0, 1)), "thrust_magn"))
    #hists.append(df.Histo1D(("thrust_costheta", "", *bins_cos_abs), "thrust_costheta"))

    df_zh = df.Filter("missingMass < 130")
    df_vbf = df.Filter("missingMass > 130")

    #########
    ### CUT 3: recoil cut (Z mass)
    #########  
    hists.append(df_zh.Histo1D(("haa_recoil_m_nOne", "", *bins_m), "haa_recoil_m"))
    if ecm == 240:
        df_zh = df_zh.Filter("haa_recoil_m > 85 && haa_recoil_m < 110")
    if ecm == 365:
        df_zh = df_zh.Filter("haa_recoil_m > 70 && haa_recoil_m < 230")
        # loose cut, as VBF has higher recoil mass
        # cut tighter later for qq/mumu channels
    hists.append(df_zh.Histo1D(("cutFlow_zh", "", *bins_count), "cut3"))


    #####
    ### CUT 4: momentum
    #####
    hists.append(df_zh.Histo1D(("haa_p_nOne", "", *bins_p), "haa_p"))
    if ecm == 240:
        df_zh = df_zh.Filter("haa_p > 20 && haa_p < 55")
    if ecm == 365:
        df_zh = df_zh.Filter("haa_p > 120 && haa_p < 150")
    hists.append(df_zh.Histo1D(("cutFlow_zh", "", *bins_count), "cut4"))


    df_zh = df_zh.Define("acoplanarity", "FCCAnalyses::acoplanarity(haa_photons)")
    df_zh = df_zh.Define("acolinearity", "FCCAnalyses::acolinearity(haa_photons)")

    ####
    ## CUT 5: acolinearity
    ####
    hists.append(df_zh.Histo1D(("acolinearity_nOne", "", *bins_aco), "acolinearity"))
    if ecm == 240:
        df_zh = df_zh.Filter("acolinearity > 0.2 && acolinearity < 0.8")
    if ecm == 365:
        ##df_zh = df_zh.Filter("acolinearity > 1.4 && acolinearity < 1.75") # very tight
        df_zh = df_zh.Filter("acolinearity > 0.8 && acolinearity < 1.75") # very tight
    hists.append(df_zh.Histo1D(("cutFlow_zh", "", *bins_count), "cut5"))


    ####
    ## CUT 6: acoplanarity
    ####
    hists.append(df_zh.Histo1D(("acoplanarity_nOne", "", *bins_aco), "acoplanarity"))
    if ecm == 240:
        df_zh = df_zh.Filter("acoplanarity > 0.05")
    if ecm == 365:
        df_zh = df_zh.Filter("acoplanarity > 0.5")
    hists.append(df_zh.Histo1D(("cutFlow_zh", "", *bins_count), "cut6"))


    ####
    ## CUT 7: cos theta(miss)
    ####
    df_zh = df_zh.Define("missingEnergy_rp", "FCCAnalyses::missingEnergy(ecm, ReconstructedParticles)")
    df_zh = df_zh.Define("missingEnergy", "missingEnergy_rp[0].energy")
    df_zh = df_zh.Define("cosThetaMiss", "FCCAnalyses::get_cosTheta_miss(missingEnergy_rp)")
    hists.append(df_zh.Histo1D(("cosThetaMiss_nOne", "", *bins_cosThetaMiss), "cosThetaMiss"))
    
    if ecm == 240:
        df_zh = df_zh.Filter("cosThetaMiss < .998")
    if ecm == 365:
        df_zh = df_zh.Filter("cosThetaMiss < .999")
    hists.append(df_zh.Histo1D(("cutFlow_zh", "", *bins_count), "cut7"))


    #########
    ### CUT 8 :cut on Higgs mass
    #########
    hists.append(df_zh.Histo1D(("haa_m_nOne", "", *bins_maa), "haa_m"))
    if ecm == 240:
        df_zh = df_zh.Filter("haa_m > 120 && haa_m < 130")
    if ecm == 365:
        df_zh = df_zh.Filter("haa_m > 120 && haa_m < 130")
    hists.append(df_zh.Histo1D(("cutFlow_zh", "", *bins_count), "cut8"))
    

    #########
    ### extra variables
    #########
    
    # missing energy
    hists.append(df_zh.Histo1D(("missingEnergy", "", *bins_m), "missingEnergy"))


    # cos theta of diphoton system
    df_zh = df_zh.Define("haa_costheta", "abs(cos(haa_theta))")
    hists.append(df_zh.Histo1D(("haa_costheta", "", *(100, 0, 1)), "haa_costheta"))
    
    # opening angle between photons
    df_zh = df_zh.Define("photon1_tlv", "ROOT::Math::PxPyPzEVector(haa_photons[0].momentum.x, haa_photons[0].momentum.y, haa_photons[0].momentum.z, haa_photons[0].energy)")
    df_zh = df_zh.Define("photon2_tlv", "ROOT::Math::PxPyPzEVector(haa_photons[1].momentum.x, haa_photons[1].momentum.y, haa_photons[1].momentum.z, haa_photons[1].energy)")
    df_zh = df_zh.Define("haa_tlv", "ROOT::Math::PxPyPzEVector(haa_cand[0].momentum.x, haa_cand[0].momentum.y, haa_cand[0].momentum.z, haa_cand[0].energy)")
    df_zh = df_zh.Define("vec3_photon1", "ROOT::Math::XYZVector(photon1_tlv.Px(),photon1_tlv.Py(),photon1_tlv.Pz());")
    df_zh = df_zh.Define("vec3_photon2", "ROOT::Math::XYZVector(photon2_tlv.Px(),photon2_tlv.Py(),photon2_tlv.Pz());")
    df_zh = df_zh.Define("vec3_diphoton", "ROOT::Math::XYZVector(haa_tlv.Px(),haa_tlv.Py(),haa_tlv.Pz());")
    df_zh = df_zh.Define("cosThetaPhotons", "abs(cos((vec3_photon1 - vec3_photon2).Theta()))")
    hists.append(df_zh.Histo1D(("cosThetaPhotons", "", *(100, 0, 1)), "cosThetaPhotons"))

    # photon energy ratio ratio
    df_zh = df_zh.Define("photon_momentum_ratio", "photon_leading_p/photon_subleading_p")
    hists.append(df_zh.Histo1D(("photon_momentum_ratio", "", *(100, 0, 10)), "photon_momentum_ratio"))

    # dphi between photons
    df_zh = df_zh.Define("photon_dphi", "ROOT::Math::VectorUtil::DeltaPhi(photon1_tlv, photon2_tlv)")
    hists.append(df_zh.Histo1D(("photon_dphi", "", *(200, -10, 10)), "photon_dphi"))




    ############################################################################


    # muons
    df_zh = df_zh.Alias("Muon0", "Muon#0.index")
    df_zh = df_zh.Define("muons_all", "FCCAnalyses::ReconstructedParticle::get(Muon0, ReconstructedParticles)")
    df_zh = df_zh.Define("muons_all_p", "FCCAnalyses::ReconstructedParticle::get_p(muons_all)")
    df_zh = df_zh.Define("muons_all_theta", "FCCAnalyses::ReconstructedParticle::get_theta(muons_all)")
    df_zh = df_zh.Define("muons_all_phi", "FCCAnalyses::ReconstructedParticle::get_phi(muons_all)")
    df_zh = df_zh.Define("muons_all_q", "FCCAnalyses::ReconstructedParticle::get_charge(muons_all)")
    df_zh = df_zh.Define("muons_all_no", "FCCAnalyses::ReconstructedParticle::get_n(muons_all)")

    df_zh = df_zh.Define("muons", "FCCAnalyses::ReconstructedParticle::sel_p(25)(muons_all)")
    df_zh = df_zh.Define("muons_p", "FCCAnalyses::ReconstructedParticle::get_p(muons)")
    df_zh = df_zh.Define("muons_theta", "FCCAnalyses::ReconstructedParticle::get_theta(muons)")
    df_zh = df_zh.Define("muons_phi", "FCCAnalyses::ReconstructedParticle::get_phi(muons)")
    df_zh = df_zh.Define("muons_q", "FCCAnalyses::ReconstructedParticle::get_charge(muons)")
    df_zh = df_zh.Define("muons_no", "FCCAnalyses::ReconstructedParticle::get_n(muons)")

    
    # electrons
    df_zh = df_zh.Alias("Electron0", "Electron#0.index")
    df_zh = df_zh.Define("electrons_all", "FCCAnalyses::ReconstructedParticle::get(Electron0, ReconstructedParticles)")
    df_zh = df_zh.Define("electrons_all_p", "FCCAnalyses::ReconstructedParticle::get_p(electrons_all)")
    df_zh = df_zh.Define("electrons_all_theta", "FCCAnalyses::ReconstructedParticle::get_theta(electrons_all)")
    df_zh = df_zh.Define("electrons_all_phi", "FCCAnalyses::ReconstructedParticle::get_phi(electrons_all)")
    df_zh = df_zh.Define("electrons_all_q", "FCCAnalyses::ReconstructedParticle::get_charge(electrons_all)")
    df_zh = df_zh.Define("electrons_all_no", "FCCAnalyses::ReconstructedParticle::get_n(electrons_all)")

    df_zh = df_zh.Define("electrons", "FCCAnalyses::ReconstructedParticle::sel_p(25)(electrons_all)")
    df_zh = df_zh.Define("electrons_p", "FCCAnalyses::ReconstructedParticle::get_p(electrons)")
    df_zh = df_zh.Define("electrons_theta", "FCCAnalyses::ReconstructedParticle::get_theta(electrons)")
    df_zh = df_zh.Define("electrons_phi", "FCCAnalyses::ReconstructedParticle::get_phi(electrons)")
    df_zh = df_zh.Define("electrons_q", "FCCAnalyses::ReconstructedParticle::get_charge(electrons)")
    df_zh = df_zh.Define("electrons_no", "FCCAnalyses::ReconstructedParticle::get_n(electrons)")


    # lepton kinematic histograms
    hists.append(df_zh.Histo1D(("muons_all_p", "", *bins_p), "muons_all_p"))
    hists.append(df_zh.Histo1D(("muons_all_theta", "", *bins_theta), "muons_all_theta"))
    hists.append(df_zh.Histo1D(("muons_all_phi", "", *bins_phi), "muons_all_phi"))
    hists.append(df_zh.Histo1D(("muons_all_q", "", *bins_charge), "muons_all_q"))
    hists.append(df_zh.Histo1D(("muons_all_no", "", *bins_count), "muons_all_no"))

    hists.append(df_zh.Histo1D(("electrons_all_p", "", *bins_p), "electrons_all_p"))
    hists.append(df_zh.Histo1D(("electrons_all_theta", "", *bins_theta), "electrons_all_theta"))
    hists.append(df_zh.Histo1D(("electrons_all_phi", "", *bins_phi), "electrons_all_phi"))
    hists.append(df_zh.Histo1D(("electrons_all_q", "", *bins_charge), "electrons_all_q"))
    hists.append(df_zh.Histo1D(("electrons_all_no", "", *bins_count), "electrons_all_no"))



    ##### CATEGORIZATION: based on #muons, # electrons, missing energy
    select_mumu = "muons_no == 2 && electrons_no == 0 && missingMass < 15"
    select_ee = "electrons_no == 2 && muons_no == 0 && missingMass < 15"
    select_nunu = "electrons_no == 0 && muons_no == 0 && missingMass > 85"
    select_qq = "electrons_no == 0 && muons_no == 0 && missingMass < 15"
    #select_tauhtauh = "electrons_no == 0 && muons_no == 0 && missingMass > 15 && missingMass < 85"




    #######
    # qq final state 
    #######
    df_qq = df_zh.Filter(select_qq)
    hists.append(df_qq.Histo1D(("cutFlow_zh", "", *bins_count), "cut9"))


    # define PF candidates collection by removing the muons
    df_qq = df_qq.Define("rps_no_photons", "FCCAnalyses::ReconstructedParticle::remove(ReconstructedParticles, photons_orig)")
    df_qq = df_qq.Define("RP_px", "FCCAnalyses::ReconstructedParticle::get_px(rps_no_photons)")
    df_qq = df_qq.Define("RP_py", "FCCAnalyses::ReconstructedParticle::get_py(rps_no_photons)")
    df_qq = df_qq.Define("RP_pz","FCCAnalyses::ReconstructedParticle::get_pz(rps_no_photons)")
    df_qq = df_qq.Define("RP_e", "FCCAnalyses::ReconstructedParticle::get_e(rps_no_photons)")
    df_qq = df_qq.Define("RP_m", "FCCAnalyses::ReconstructedParticle::get_mass(rps_no_photons)")
    df_qq = df_qq.Define("RP_q", "FCCAnalyses::ReconstructedParticle::get_charge(rps_no_photons)")
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

    hists.append(df_qq.Histo1D(("zqq_jet1_p_nOne", "", *bins_p), "jet1_p"))
    hists.append(df_qq.Histo1D(("zqq_jet2_p_nOne", "", *bins_p), "jet2_p"))
    if ecm == 240:
        pass
        #df_qq = df_qq.Filter("dijet_p > 25 && dijet_p < 55")
    if ecm == 365:
        df_qq = df_qq.Filter("jet1_p > 20 && jet2_p > 10")



    hists.append(df_qq.Histo1D(("zqq_qq_p_nOne", "", *bins_p), "dijet_p"))
    if ecm == 240:
        df_qq = df_qq.Filter("dijet_p > 25 && dijet_p < 55")
    if ecm == 365:
        df_qq = df_qq.Filter("dijet_p > 100 && dijet_p < 150")

    hists.append(df_qq.Histo1D(("zqq_m_nOne", "", *bins_m), "dijet_m"))
    if ecm == 240:
        df_qq = df_qq.Filter("dijet_m < 105 && dijet_m > 75")
    if ecm == 365:
        df_qq = df_qq.Filter("dijet_m < 130 && dijet_m > 80")
    

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
    df_nunu = df_zh.Filter(select_nunu)
    hists.append(df_nunu.Histo1D(("cutFlow_zh", "", *bins_count), "cut10"))
    hists.append(df_nunu.Histo1D(("znunu_haa_m", "", *bins_m_zoom), "haa_m"))

    # extra variables
    hists.append(df_nunu.Histo1D(("znunu_cosThetaMiss", "", *bins_cosThetaMiss), "cosThetaMiss"))
    hists.append(df_nunu.Histo1D(("znunu_missingEnergy", "", *bins_m), "missingEnergy"))
    hists.append(df_nunu.Histo1D(("znunu_cosThetaPhotons", "", *(100, 0, 1)), "cosThetaPhotons"))



    #######
    # mumu final state
    #######
    
    ## can use the recoil as observable? Better sensitivity

    df_mumu = df_zh.Filter(select_mumu)
    hists.append(df_mumu.Histo1D(("cutFlow_zh", "", *bins_count), "cut11"))

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


    # construct recoil around tight aa mass
    # recoil has better resolution than m(aa)
    #df_mumu = df_mumu.Filter("haa_m > 123 && haa_m < 127")
    #df_mumu = df_mumu.Define("init_tlv", f"TLorentzVector ret; ret.SetPxPyPzE(0, 0, 0, {ecm}); return ret;")
    #df_mumu = df_mumu.Define("zmumu_recoil_tlv", "init_tlv-dimuon")
    #df_mumu = df_mumu.Define("zmumu_recoil_m", "zmumu_recoil_tlv.M()")
    #hists.append(df_mumu.Histo1D(("zmumu_recoil_m", "", *(100, 120, 130)), "zmumu_recoil_m"))



    #######
    # ee final state
    #######
    
    ## can use the recoil as observable? Better sensitivity
    df_ee = df_zh.Filter(select_ee)
    hists.append(df_ee.Histo1D(("cutFlow_zh", "", *bins_count), "cut12"))

    df_ee = df_ee.Define("electrons_tlv", "FCCAnalyses::makeLorentzVectors(electrons)")
    df_ee = df_ee.Define("dielectron", "electrons_tlv[0]+electrons_tlv[1]")
    df_ee = df_ee.Define("dielectron_m", "dielectron.M()")

    hists.append(df_ee.Histo1D(("zee_m_nOne", "", *bins_m), "dielectron_m"))
    df_ee = df_ee.Filter("dielectron_m > 80 && dielectron_m < 100")
    hists.append(df_ee.Histo1D(("zee_m", "", *bins_m), "dielectron_m"))
    hists.append(df_ee.Histo1D(("zee_haa_m", "", *bins_m_zoom), "haa_m"))






    ### VBF ANALYSIS

    #########
    ### CUT 3: recoil cut (Z mass)
    #########
    hists.append(df_vbf.Histo1D(("vbf_haa_recoil_m_nOne", "", *bins_m), "haa_recoil_m"))
    df_vbf = df_vbf.Filter("haa_recoil_m < 240")
    hists.append(df_vbf.Histo1D(("cutFlow_vbf", "", *bins_count), "cut3"))

    #####
    ### CUT 4: momentum
    #####
    
    hists.append(df_vbf.Histo1D(("vbf_haa_p_nOne", "", *bins_p), "haa_p"))
    ##df_vbf = df_vbf.Filter("haa_p > 20")
    df_vbf = df_vbf.Filter("haa_p > 30 && haa_p < 130")
    hists.append(df_vbf.Histo1D(("cutFlow_vbf", "", *bins_count), "cut4"))


    df_vbf = df_vbf.Define("acoplanarity", "FCCAnalyses::acoplanarity(haa_photons)")
    df_vbf = df_vbf.Define("acolinearity", "FCCAnalyses::acolinearity(haa_photons)")

    ####
    ## CUT 5: acolinearity
    ####
    hists.append(df_vbf.Histo1D(("vbf_acolinearity_nOne", "", *bins_aco), "acolinearity"))
    #df_vbf = df_vbf.Filter("acolinearity > 0.2 && acolinearity < 1.6")
    df_vbf = df_vbf.Filter("acolinearity > 0.35 && acolinearity < 1.6")
    hists.append(df_vbf.Histo1D(("cutFlow_vbf", "", *bins_count), "cut5"))


    ####
    ## CUT 6: acoplanarity
    ####
    hists.append(df_vbf.Histo1D(("vbf_acoplanarity_nOne", "", *bins_aco), "acoplanarity"))
    #df_vbf = df_vbf.Filter("acoplanarity > 0.1")
    df_vbf = df_vbf.Filter("acoplanarity > 0.2")
    hists.append(df_vbf.Histo1D(("cutFlow_vbf", "", *bins_count), "cut6"))

    ####
    ## CUT 7: missing mass
    ####
    hists.append(df_vbf.Histo1D(("vbf_missingMass_nOne", "", *bins_m), "missingMass"))
    df_vbf = df_vbf.Filter("missingMass > 170")
    hists.append(df_vbf.Histo1D(("cutFlow_vbf", "", *bins_count), "cut7"))


    #########
    ### CUT 8 :cut on Higgs mass
    #########
    hists.append(df_vbf.Histo1D(("vbf_haa_m_nOne", "", *bins_maa), "haa_m"))
    df_vbf = df_vbf.Filter("haa_m > 120 && haa_m < 130")
    hists.append(df_vbf.Histo1D(("cutFlow_vbf", "", *bins_count), "cut8"))

    hists.append(df_vbf.Histo1D(("vbf_haa_m", "", *bins_m_zoom), "haa_m"))


    df_vbf = df_vbf.Define("missingEnergy_rp", "FCCAnalyses::missingEnergy(ecm, ReconstructedParticles)")
    df_vbf = df_vbf.Define("missingEnergy", "missingEnergy_rp[0].energy")
    df_vbf = df_vbf.Define("cosThetaMiss", "FCCAnalyses::get_cosTheta_miss(missingEnergy_rp)")
    hists.append(df_vbf.Histo1D(("vbf_cosThetaMiss_nOne", "", *bins_cosThetaMiss), "cosThetaMiss"))


    # cos theta of diphoton system
    df_vbf = df_vbf.Define("haa_costheta", "abs(cos(haa_theta))")
    hists.append(df_vbf.Histo1D(("vbf_haa_costheta_nOne", "", *(100, 0, 1)), "haa_costheta"))

    return hists, weightsum
