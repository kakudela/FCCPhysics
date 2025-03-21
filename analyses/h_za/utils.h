

namespace FCCAnalyses {


Vec_i jetTruthFinder(std::vector<std::vector<int>> constituents, Vec_rp reco, Vec_mc mc, Vec_i mcind) {
    // jet truth=finder: match the gen-level partons (eventually with gluons) with the jet constituents
    // matching by mimimizing the sum of dr of the parton and all the jet constituents 

    Vec_tlv genQuarks; // Lorentz-vector of potential partons (gen-level)
    Vec_i genQuarks_pdgId; // corresponding PDG ID
    for(size_t i = 0; i < mc.size(); ++i) {
        int pdgid = abs(mc.at(i).PDG);
        if(pdgid > 6) continue; // only quarks 
        //if(pdgid > 6 and pdgid != 21) continue; // only quarks and gluons
        TLorentzVector tlv;
        tlv.SetXYZM(mc.at(i).momentum.x,mc.at(i).momentum.y,mc.at(i).momentum.z,mc.at(i).mass);
        genQuarks.push_back(tlv);
        genQuarks_pdgId.push_back(mc.at(i).PDG);
    }

    Vec_tlv recoParticles; // Lorentz-vector of all reconstructed particles
    for(size_t i = 0; i < reco.size(); ++i) {
        auto & p = reco[i];
        TLorentzVector tlv;
        tlv.SetXYZM(p.momentum.x, p.momentum.y, p.momentum.z, p.mass);
        recoParticles.push_back(tlv);
    }

    Vec_i usedIdx;
    Vec_i result;
    for(size_t iJet = 0; iJet < constituents.size(); ++iJet) {
        Vec_d dr;
        for(size_t iGen = 0; iGen < genQuarks.size(); ++iGen) {
            if(std::find(usedIdx.begin(), usedIdx.end(), iGen) != usedIdx.end()) {
                dr.push_back(1e99); // set infinite dr, skip
                continue;
            }
            dr.push_back(0);
            for(size_t i = 0; i < constituents[iJet].size(); ++i) {
                dr[iGen] += recoParticles[constituents[iJet][i]].DeltaR(genQuarks[iGen]);
            }
        }
        int maxDrIdx = std::min_element(dr.begin(),dr.end()) - dr.begin();
        usedIdx.push_back(maxDrIdx);
        result.push_back(genQuarks_pdgId[maxDrIdx]);

    }
    return result;
}

float pairing_H_Z_a_chi2(float ecm, TLorentzVector Z1, TLorentzVector Z2, TLorentzVector a) {

    Vec_tlv ret; // returns Higgs and associated Z

    TLorentzVector init;
    init.SetPxPyPzE(0, 0, 0, ecm);

    auto Z1a = Z1 + a;
    auto Z2a = Z2 + a;

    auto Z1_rec = init - Z1;
    auto Z2_rec = init - Z2;
    
    //float chi2 = std::pow(qqa_m-125, 2) + 0.0*std::pow(vv_recoil_m-125, 2) - std::pow(qq_recoil_m-125, 2) - std::pow(vva_m-125, 2);
    float chi2 = std::pow(Z1a.M()-125, 2) + std::pow(Z2_rec.M()-125, 2) - std::pow(Z1_rec.M()-125, 2) - std::pow(Z2a.M()-125, 2);

    return chi2;
}


Vec_tlv pairing_H_Z_a(float ecm, TLorentzVector Z1, TLorentzVector Z2, TLorentzVector a) {

    Vec_tlv ret; // returns Higgs and associated Z

    TLorentzVector init;
    init.SetPxPyPzE(0, 0, 0, ecm);

    auto Z1a = Z1 + a;
    auto Z2a = Z2 + a;

    auto Z1_rec = init - Z1;
    auto Z2_rec = init - Z2;


    // case H=Z1a
    float chi2_1 = std::pow(Z2_rec.M()-125, 2) + std::pow(Z1a.M()-125, 2) + std::pow(Z1a.P()-50, 2);

    // case H=Z2a
    float chi2_2 = std::pow(Z1_rec.M()-125, 2) + std::pow(Z2a.M()-125, 2) + std::pow(Z2a.P()-50, 2);

    if(chi2_1 < chi2_2) {
        ret.push_back(Z1a); // Higgs
        ret.push_back(Z2); // associated Z
        ret.push_back(Z1); // Z from Higgs
    }
    else {
        ret.push_back(Z2a); // Higgs
        ret.push_back(Z1); // associated Z
        ret.push_back(Z2); // Z from Higgs
    }

    return ret;
}




Vec_tlv pairing_Z(Vec_tlv in) {
    // assume 4 input jets
    Vec_tlv ret;

    TLorentzVector Z1, Z2;
    auto j1 = in[0];
    auto j2 = in[1];
    auto j3 = in[2];
    auto j4 = in[3];

    float chi2_1 = std::pow((j1+j2).M()-91.2, 2) + std::pow((j3+j4).M()-91.2, 2);
    float chi2_2 = std::pow((j1+j3).M()-91.2, 2) + std::pow((j2+j4).M()-91.2, 2);
    float chi2_3 = std::pow((j1+j4).M()-91.2, 2) + std::pow((j2+j3).M()-91.2, 2);

    if(chi2_1<chi2_2 && chi2_1<chi2_3) {
        Z1 = j1+j2;
        Z2 = j3+j4;
    }
    else if(chi2_2<chi2_1 && chi2_2<chi2_3) {
        Z1 = j1+j3;
        Z2 = j2+j4;
    }
    else if(chi2_3<chi2_1 && chi2_3<chi2_2) {
        Z1 = j1+j4;
        Z2 = j2+j3;
    }

    ret.push_back(Z1);
    ret.push_back(Z2);
    return ret;

}




}