#pragma once
#include <cmath>
#include <algorithm>
#include <vector>
#include <numeric>
#include <ROOT/RVec.hxx>
#include <unordered_map>
#include "edm4hep/SimTrackerHitData.h"
#include "edm4hep/MCParticleData.h"

using Vec_b  = ROOT::VecOps::RVec<bool>;
using Vec_f  = ROOT::VecOps::RVec<float>;
using Vec_d  = ROOT::VecOps::RVec<double>;
using Vec_i  = ROOT::VecOps::RVec<int>;
using Vec_ui = ROOT::VecOps::RVec<unsigned int>;

using Vec_SimTrackerHitData = ROOT::VecOps::RVec<edm4hep::SimTrackerHitData>;
using Vec_MCParticleData    = ROOT::VecOps::RVec<edm4hep::MCParticleData>;

// sim hit helpers
Vec_f getSimHitPosition_x(const Vec_SimTrackerHitData& in){
    Vec_f ret; ret.reserve(in.size());
    for (const auto& h : in) ret.push_back(h.position[0]);
    return ret;
}
Vec_f getSimHitPosition_y(const Vec_SimTrackerHitData& in){
    Vec_f ret; ret.reserve(in.size());
    for (const auto& h : in) ret.push_back(h.position[1]);
    return ret;
}
Vec_f getSimHitPosition_z(const Vec_SimTrackerHitData& in){
    Vec_f ret; ret.reserve(in.size());
    for (const auto& h : in) ret.push_back(h.position[2]);
    return ret;
}
Vec_f getSimHitPosition_r(const Vec_SimTrackerHitData& in){
    Vec_f ret; ret.reserve(in.size());
    for (const auto& h : in){
        const float x=h.position[0], y=h.position[1];
        ret.push_back(std::sqrt(x*x + y*y));
    }
    return ret;
}
Vec_f getSimHitPosition_theta(const Vec_SimTrackerHitData& in, bool deg=false){
    Vec_f ret; ret.reserve(in.size());
    for (const auto& h : in){
        const float x=h.position[0], y=h.position[1], z=h.position[2];
        const double r3 = std::sqrt(double(x)*x + double(y)*y + double(z)*z);
        double th = (r3>0.0) ? std::acos(z / r3) : 0.0;
        if (deg) th *= 180.0 / M_PI;
        ret.push_back((float)th);
    }
    return ret;
}
Vec_f getSimHitPosition_phi(const Vec_SimTrackerHitData& in, bool deg=false){
    Vec_f ret; ret.reserve(in.size());
    for (const auto& h : in){
        double ph = std::atan2(h.position[1], h.position[0]);
        if (deg) ph *= 180.0 / M_PI;
        ret.push_back((float)ph);
    }
    return ret;
}
Vec_b isProducedBySecondary(const Vec_SimTrackerHitData& in){
    Vec_b ret; ret.reserve(in.size());
    for (const auto& h : in) ret.push_back( (h.quality & (1<<30)) != 0 );
    return ret;
}
Vec_f getEnergyDeposition(const Vec_SimTrackerHitData& in){
    Vec_f ret; ret.reserve(in.size());
    for (const auto& h : in) ret.push_back(h.eDep);
    return ret;
}
Vec_i getCellID(const Vec_SimTrackerHitData& in){
    Vec_i ret; ret.reserve(in.size());
    for (const auto& h : in) ret.push_back((int)h.cellID);
    return ret;
}
Vec_i getSimHitLayer(const Vec_f& hits_r, const Vec_f& radii, float tol_mm=4.f){
    Vec_i ret; ret.reserve(hits_r.size());
    for (float r : hits_r){
        int layer = -1;
        for (int i=0;i<(int)radii.size();++i){
            if (std::abs(r - radii[i]) < tol_mm){ layer = i; break; }
        }
        ret.push_back(layer);
    }
    return ret;
}

// mc helpers
Vec_MCParticleData getMCParticle(const Vec_SimTrackerHitData& in, const Vec_MCParticleData& mc, const Vec_i& idx){
    Vec_MCParticleData ret; ret.reserve(in.size());
    for (int i=0;i<(int)in.size();++i){
        if (i < (int)idx.size() && idx[i] >= 0 && idx[i] < (int)mc.size()) ret.push_back(mc[idx[i]]);
        else ret.push_back(edm4hep::MCParticleData());
    }
    return ret;
}
Vec_f getMCVertex_x(const Vec_MCParticleData& in){
    Vec_f ret; ret.reserve(in.size());
    for (const auto& p : in) ret.push_back((float)p.vertex[0]);
    return ret;
}
Vec_f getMCVertex_y(const Vec_MCParticleData& in){
    Vec_f ret; ret.reserve(in.size());
    for (const auto& p : in) ret.push_back((float)p.vertex[1]);
    return ret;
}
Vec_f getMCVertex_z(const Vec_MCParticleData& in){
    Vec_f ret; ret.reserve(in.size());
    for (const auto& p : in) ret.push_back((float)p.vertex[2]);
    return ret;
}
Vec_f getMCVertex_r(const Vec_MCParticleData& in){
    Vec_f ret; ret.reserve(in.size());
    for (const auto& p : in){
        const double x=p.vertex[0], y=p.vertex[1];
        ret.push_back((float)std::sqrt(x*x + y*y));
    }
    return ret;
}
Vec_f getMCMomentum_px(const Vec_MCParticleData& in){
    Vec_f ret; ret.reserve(in.size());
    for (const auto& p : in) ret.push_back((float)p.momentum[0]);
    return ret;
}
Vec_f getMCMomentum_py(const Vec_MCParticleData& in){
    Vec_f ret; ret.reserve(in.size());
    for (const auto& p : in) ret.push_back((float)p.momentum[1]);
    return ret;
}
Vec_f getMCMomentum_pz(const Vec_MCParticleData& in){
    Vec_f ret; ret.reserve(in.size());
    for (const auto& p : in) ret.push_back((float)p.momentum[2]);
    return ret;
}
Vec_f getMCMomentum_theta(const Vec_MCParticleData& in, bool deg=false){
    Vec_f ret; ret.reserve(in.size());
    for (const auto& p : in){
        const double px=p.momentum[0], py=p.momentum[1], pz=p.momentum[2];
        const double p3 = std::sqrt(px*px + py*py + pz*pz);
        double th = (p3>0.0) ? std::acos(pz / p3) : 0.0;
        if (deg) th *= 180.0 / M_PI;
        ret.push_back((float)th);
    }
    return ret;
}
Vec_f getMCMomentum_phi(const Vec_MCParticleData& in, bool deg=false){
    Vec_f ret; ret.reserve(in.size());
    for (const auto& p : in){
        double ph = std::atan2(p.momentum[1], p.momentum[0]);
        if (deg) ph *= 180.0 / M_PI;
        ret.push_back((float)ph);
    }
    return ret;
}
Vec_f getMCMomentum_pt(const Vec_MCParticleData& in){
    Vec_f ret; ret.reserve(in.size());
    for (const auto& p : in){
        const double px=p.momentum[0], py=p.momentum[1];
        ret.push_back((float)std::sqrt(px*px + py*py));
    }
    return ret;
}
Vec_f getMCMomentum_p(const Vec_MCParticleData& in){
    Vec_f ret; ret.reserve(in.size());
    for (const auto& p : in){
        const double px=p.momentum[0], py=p.momentum[1], pz=p.momentum[2];
        ret.push_back((float)std::sqrt(px*px + py*py + pz*pz));
    }
    return ret;
}
Vec_f getMCEnergy(const Vec_MCParticleData& in){
    Vec_f ret; ret.reserve(in.size());
    for (const auto& p : in){
        const double px=p.momentum[0], py=p.momentum[1], pz=p.momentum[2];
        const double p3 = std::sqrt(px*px + py*py + pz*pz);
        const double m  = (double)p.mass;
        const double E  = std::sqrt(p3*p3 + m*m);
        ret.push_back((float)E);
    }
    return ret;
}

Vec_i getMCGeneratorStatus(const Vec_MCParticleData& in){
    Vec_i ret; ret.reserve(in.size());
    for (const auto& p : in) ret.push_back((int)p.generatorStatus);
    return ret;
}
Vec_i getMCPDGID(const Vec_MCParticleData& in, bool absval=false){
    Vec_i ret; ret.reserve(in.size());
    for (const auto& p : in){
        int id = (int)p.PDG;
        ret.push_back(absval ? std::abs(id) : id);
    }
    return ret;
}

// vector math helpers
Vec_f cosFromRadians(const Vec_f& th){
    Vec_f out; out.reserve(th.size());
    for(auto x : th) out.push_back((float)std::cos(x));
    return out;
}
Vec_f vecAbs(const Vec_f& v){
    Vec_f out; out.reserve(v.size());
    for(auto x : v) out.push_back(std::fabs(x));
    return out;
}
Vec_f log10Vec(const Vec_f& v){
    Vec_f out; out.reserve(v.size());
    for(auto x : v) out.push_back((float)std::log10(std::max(1e-12f, x)));
    return out;
}

// angle wrapping helpers (degrees)
inline float wrapDeg360(float a){
    float x = std::fmod(a, 360.0f);
    if (x < 0.0f) x += 360.0f;
    return x;
}
template <typename T>
ROOT::VecOps::RVec<float> wrapVecDeg360(const ROOT::VecOps::RVec<T>& v){
    ROOT::VecOps::RVec<float> out; out.reserve(v.size());
    for (const auto& a : v) out.push_back(wrapDeg360((float)a));
    return out;
}

// per-MC hit multiplicity helpers
Vec_i mcHitMultiplicity(const Vec_i& mcidx){
    std::unordered_map<int,int> m;
    for (int id : mcidx){
        if (id < 0) continue;
        m[id] += 1;
    }
    Vec_i out;
    out.reserve(m.size());
    for (const auto& kv : m) out.push_back(kv.second);
    return out;
}
int mcUniqueCount(const Vec_i& mcidx){
    std::unordered_map<int,int> m;
    for (int id : mcidx){
        if (id < 0) continue;
        m[id] = 1;
    }
    return (int)m.size();
}

// consecutive deltas between hits of same MC
Vec_f deltaSeqSameMC(const Vec_i& mcidx, const Vec_f& vals){
    const int n = std::min((int)mcidx.size(), (int)vals.size());
    Vec_f out;
    if (n <= 1) return out;
    out.reserve(n-1);
    for (int i=1;i<n;i++){
        if (mcidx[i] == mcidx[i-1] && mcidx[i] >= 0) out.push_back(vals[i] - vals[i-1]);
    }
    return out;
}
Vec_f deltaPhiSeqSameMC(const Vec_i& mcidx, const Vec_f& phi_deg){
    const int n = std::min((int)mcidx.size(), (int)phi_deg.size());
    Vec_f out;
    if (n <= 1) return out;
    out.reserve(n-1);
    for (int i=1;i<n;i++){
        if (mcidx[i] != mcidx[i-1] || mcidx[i] < 0) continue;
        float d = phi_deg[i] - phi_deg[i-1];
        while (d > 180.0f) d -= 360.0f;
        while (d < -180.0f) d += 360.0f;
        out.push_back(d);
    }
    return out;
}

// SimHit momentum, MC beta, de-dup, incidence angle helpers

Vec_f getSimHitMomentum_x(const Vec_SimTrackerHitData& in){
    Vec_f ret; ret.reserve(in.size());
    for (const auto& h : in) ret.push_back((float)h.momentum[0]);
    return ret;
}
Vec_f getSimHitMomentum_y(const Vec_SimTrackerHitData& in){
    Vec_f ret; ret.reserve(in.size());
    for (const auto& h : in) ret.push_back((float)h.momentum[1]);
    return ret;
}
Vec_f getSimHitMomentum_z(const Vec_SimTrackerHitData& in){
    Vec_f ret; ret.reserve(in.size());
    for (const auto& h : in) ret.push_back((float)h.momentum[2]);
    return ret;
}
Vec_f getSimHitMomentum_pt(const Vec_SimTrackerHitData& in){
    Vec_f ret; ret.reserve(in.size());
    for (const auto& h : in){
        const float px=(float)h.momentum[0], py=(float)h.momentum[1];
        ret.push_back(std::sqrt(px*px + py*py));
    }
    return ret;
}
Vec_f getSimHitMomentum_p(const Vec_SimTrackerHitData& in){
    Vec_f ret; ret.reserve(in.size());
    for (const auto& h : in){
        const float px=(float)h.momentum[0], py=(float)h.momentum[1], pz=(float)h.momentum[2];
        ret.push_back(std::sqrt(px*px + py*py + pz*pz));
    }
    return ret;
}

Vec_f getMCBeta(const Vec_MCParticleData& in){
    Vec_f ret; ret.reserve(in.size());
    for (const auto& p : in){
        const double px=p.momentum[0], py=p.momentum[1], pz=p.momentum[2];
        const double p3 = std::sqrt(px*px + py*py + pz*pz);
        const double m  = (double)p.mass;
        const double E  = std::sqrt(p3*p3 + m*m);
        ret.push_back((float)((E>0.0) ? (p3 / E) : 0.0));
    }
    return ret;
}


int mergedHitCountConsecutiveByMC(const Vec_i& mcidx, const Vec_f& x, const Vec_f& y, const Vec_f& z, float radius_mm){
    const float r2 = radius_mm * radius_mm;
    int kept = 0;
    std::unordered_map<int, int> has_prev;
    std::unordered_map<int, float> px;
    std::unordered_map<int, float> py;
    std::unordered_map<int, float> pz;

    for(int i=0;i<(int)mcidx.size();++i){
        const int id = mcidx[i];
        if(id < 0) continue;

        auto it = has_prev.find(id);
        if(it == has_prev.end()){
            has_prev[id] = 1;
            px[id] = x[i]; py[id] = y[i]; pz[id] = z[i];
            kept += 1;
            continue;
        }

        const float dx = x[i] - px[id];
        const float dy = y[i] - py[id];
        const float dz = z[i] - pz[id];
        const float d2 = dx*dx + dy*dy + dz*dz;
        if(d2 < r2){
            continue;
        }

        px[id] = x[i]; py[id] = y[i]; pz[id] = z[i];
        kept += 1;
    }
    return kept;
}

Vec_i mcHitMultiplicityMergedConsecutiveByMC(const Vec_i& mcidx, const Vec_f& x, const Vec_f& y, const Vec_f& z, float radius_mm){
    const float r2 = radius_mm * radius_mm;

    std::unordered_map<int, int> count;
    std::unordered_map<int, int> has_prev;
    std::unordered_map<int, float> px;
    std::unordered_map<int, float> py;
    std::unordered_map<int, float> pz;

    for(int i=0;i<(int)mcidx.size();++i){
        const int id = mcidx[i];
        if(id < 0) continue;

        auto it = has_prev.find(id);
        if(it == has_prev.end()){
            has_prev[id] = 1;
            px[id] = x[i]; py[id] = y[i]; pz[id] = z[i];
            count[id] += 1;
            continue;
        }

        const float dx = x[i] - px[id];
        const float dy = y[i] - py[id];
        const float dz = z[i] - pz[id];
        const float d2 = dx*dx + dy*dy + dz*dz;
        if(d2 < r2){
            continue;
        }

        px[id] = x[i]; py[id] = y[i]; pz[id] = z[i];
        count[id] += 1;
    }

    Vec_i out;
    out.reserve(count.size());
    for(const auto& kv : count) out.push_back(kv.second);
    return out;
}

Vec_f getBarrelIncidenceAngleDeg(const Vec_f& x, const Vec_f& y, const Vec_f& px_in, const Vec_f& py_in, const Vec_f& pz_in){
    const int n = std::min({(int)x.size(), (int)y.size(), (int)px_in.size(), (int)py_in.size(), (int)pz_in.size()});
    Vec_f out; out.reserve(n);

    for(int i=0;i<n;++i){
        const double nx = (double)x[i];
        const double ny = (double)y[i];
        const double nr = std::sqrt(nx*nx + ny*ny);

        const double px = (double)px_in[i];
        const double py = (double)py_in[i];
        const double pz = (double)pz_in[i];
        const double pm = std::sqrt(px*px + py*py + pz*pz);

        if(nr <= 0.0 || pm <= 0.0){
            out.push_back(0.0f);
            continue;
        }

        const double nhx = nx / nr;
        const double nhy = ny / nr;

        const double phx = px / pm;
        const double phy = py / pm;

        double c = nhx*phx + nhy*phy;
        c = std::max(-1.0, std::min(1.0, c));
        const double ang = std::acos(c) * 180.0 / M_PI;
        out.push_back((float)ang);
    }
    return out;
}

Vec_f getBarrelIncidenceAnglePhiDeg(const Vec_f& x, const Vec_f& y, const Vec_f& px_in, const Vec_f& py_in, const Vec_f& pz_in){
    const int n = std::min({(int)x.size(), (int)y.size(), (int)px_in.size(), (int)py_in.size(), (int)pz_in.size()});
    Vec_f out; out.reserve(n);

    for(int i=0;i<n;++i){
        const double nx = (double)x[i];
        const double ny = (double)y[i];
        const double nr = std::sqrt(nx*nx + ny*ny);

        const double px = (double)px_in[i];
        const double py = (double)py_in[i];
        const double pt = std::sqrt(px*px + py*py);

        if(nr <= 0.0 || pt <= 0.0){
            out.push_back(0.0f);
            continue;
        }

        const double nhx = nx / nr;
        const double nhy = ny / nr;

        const double phx = px / pt;
        const double phy = py / pt;

        double c = nhx*phx + nhy*phy;
        c = std::max(-1.0, std::min(1.0, c));
        const double ang = std::acos(c) * 180.0 / M_PI;
        out.push_back((float)ang);
    }
    return out;
}

Vec_f getBarrelIncidenceAngleThetaDeg(const Vec_f& x, const Vec_f& y, const Vec_f& px_in, const Vec_f& py_in, const Vec_f& pz_in){
    const int n = std::min({(int)x.size(), (int)y.size(), (int)px_in.size(), (int)py_in.size(), (int)pz_in.size()});
    Vec_f out; out.reserve(n);

    for(int i=0;i<n;++i){
        const double nx = (double)x[i];
        const double ny = (double)y[i];
        const double nr = std::sqrt(nx*nx + ny*ny);

        const double px = (double)px_in[i];
        const double py = (double)py_in[i];
        const double pz = (double)pz_in[i];

        if(nr <= 0.0){
            out.push_back(0.0f);
            continue;
        }

        // radial component of momentum (projection onto r-hat)
        const double pr = (nx*px + ny*py) / nr;

        // momentum projected into r-z plane
        const double prm = std::sqrt(pr*pr + pz*pz);
        if(prm <= 0.0){
            out.push_back(0.0f);
            continue;
        }

        // angle between +r direction and (pr, pz)
        double c = pr / prm;
        c = std::max(-1.0, std::min(1.0, c));
        const double ang = std::acos(c) * 180.0 / M_PI;
        out.push_back((float)ang);
    }
    return out;
}


// Elementwise safe division: out[i] = (den[i] != 0) ? num[i]/den[i] : 0
Vec_f vecDiv(const Vec_f& num, const Vec_f& den) {
    Vec_f out;
    const int n = (int)num.size();
    if ((int)den.size() != n) return out;
    out.reserve(n);
    for (int i = 0; i < n; ++i) {
        const float d = den[i];
        out.push_back((std::abs(d) > 0.0f) ? (num[i] / d) : 0.0f);
    }
    return out;
}

Vec_b maskMCSelectedHasHit(const Vec_ui& mc_sel_idx, const Vec_i& hit_mcidx){
    std::unordered_map<unsigned int, int> has;
    for (int id : hit_mcidx){
        if (id >= 0) has[(unsigned int)id] = 1;
    }
    Vec_b out;
    out.reserve(mc_sel_idx.size());
    for (unsigned int idx : mc_sel_idx){
        out.push_back(has.find(idx) != has.end());
    }
    return out;
}
