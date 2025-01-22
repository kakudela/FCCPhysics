

namespace FCCAnalyses {

Vec_i ww_decay_mode(Vec_mc mc, Vec_i ind) {
   Vec_i res; // returns vector of 4 indices (first two PDG of W+ daughters, second two PDG of W-)
   res.push_back(-99);
   res.push_back(-99);
   res.push_back(-99);
   res.push_back(-99);
   for(size_t i = 0; i < mc.size(); ++i) {
        auto & p = mc[i];
        if(std::abs(p.PDG) != 24) continue;

        int ds = p.daughters_begin;
        int de = p.daughters_end;
        int idx_ds = ind[ds];
        int idx_de = ind[de-1];
        int pdg_d1 = mc[idx_ds].PDG;
        int pdg_d2 = mc[idx_de].PDG;

        if(std::abs(pdg_d1) == 24 or std::abs(pdg_d2) == 24) continue;
        if(p.PDG == 24) {
            res[0] = pdg_d1;
            res[1] = pdg_d2;
        }
        else {
            res[2] = pdg_d1;
            res[3] = pdg_d2;
        }
   }
   return res;
}

}