#pragma once
#include <cmath>
#include "edm4hep/SimTrackerHitData.h"
#include "edm4hep/MCParticleData.h"
#include "ROOT/RVec.hxx"

using Vec_b  = ROOT::VecOps::RVec<bool>;
using Vec_i  = ROOT::VecOps::RVec<int>;
using Vec_ui = ROOT::VecOps::RVec<unsigned int>;
using Vec_f  = ROOT::VecOps::RVec<float>;
using Vec_d  = ROOT::VecOps::RVec<double>;

using Vec_SimHit = ROOT::VecOps::RVec<edm4hep::SimTrackerHitData>;
using Vec_MCPart = ROOT::VecOps::RVec<edm4hep::MCParticleData>;

/* ------------------------- Sim-tracker-hits ------------------------- */
inline Vec_f getSimHitPosition_x(const Vec_SimHit& in){ Vec_f o(in.size()); for(size_t i=0;i<in.size();++i) o[i]=in[i].position[0]; return o; }
inline Vec_f getSimHitPosition_y(const Vec_SimHit& in){ Vec_f o(in.size()); for(size_t i=0;i<in.size();++i) o[i]=in[i].position[1]; return o; }
inline Vec_f getSimHitPosition_z(const Vec_SimHit& in){ Vec_f o(in.size()); for(size_t i=0;i<in.size();++i) o[i]=in[i].position[2]; return o; }
inline Vec_f getSimHitPosition_r(const Vec_SimHit& in){ Vec_f o(in.size()); for(size_t i=0;i<in.size();++i) o[i]=std::hypot(in[i].position[0],in[i].position[1]); return o; }

inline Vec_f getSimHitPosition_theta(const Vec_SimHit& in,bool deg=false){
    Vec_f o(in.size());
    for(size_t i=0;i<in.size();++i){
        const auto& h=in[i];
        double r=std::sqrt(h.position[0]*h.position[0]+h.position[1]*h.position[1]+h.position[2]*h.position[2]);
        double t=r>0?std::acos(h.position[2]/r):0.; if(deg) t*=180./M_PI; o[i]=t;
    } return o;
}
inline Vec_f getSimHitPosition_phi(const Vec_SimHit& in,bool deg=false){
    Vec_f o(in.size());
    for(size_t i=0;i<in.size();++i){ double p=std::atan2(in[i].position[1],in[i].position[0]); if(deg) p*=180./M_PI; o[i]=p; }
    return o;
}

inline Vec_f getSimHitEDep(const Vec_SimHit& in){ Vec_f o(in.size()); for(size_t i=0;i<in.size();++i) o[i]=in[i].EDep*1e6f; return o; }   // GeV→keV
inline Vec_f getSimHitPathLength(const Vec_SimHit& in){ Vec_f o(in.size()); for(size_t i=0;i<in.size();++i) o[i]=in[i].pathLength; return o; }
inline Vec_ui getSimHitCellID(const Vec_SimHit& in){ Vec_ui o(in.size()); for(size_t i=0;i<in.size();++i) o[i]=(unsigned)in[i].cellID; return o; }

inline Vec_b filterSimHitsPrimary(const Vec_SimHit& in){
    Vec_b o(in.size());
    for(size_t i=0;i<in.size();++i){ bool secondary=(in[i].quality>>30)&1; o[i]=!secondary; }
    return o;
}
inline Vec_i getSimHitLayer(const Vec_f& r,const Vec_f& radii){
    Vec_i o(r.size());
    for(size_t i=0;i<r.size();++i){ int L=-1; for(size_t k=0;k<radii.size();++k) if(std::fabs(r[i]-radii[k])<4.f){L=k;break;} o[i]=L; }
    return o;
}

/* ------------------------------ MC-truth ----------------------------- */
inline Vec_d getMCParticleVertex_x(const Vec_MCPart& in){ Vec_d o(in.size()); for(size_t i=0;i<in.size();++i) o[i]=in[i].vertex[0]; return o; }
inline Vec_d getMCParticleVertex_y(const Vec_MCPart& in){ Vec_d o(in.size()); for(size_t i=0;i<in.size();++i) o[i]=in[i].vertex[1]; return o; }
inline Vec_d getMCParticleVertex_z(const Vec_MCPart& in){ Vec_d o(in.size()); for(size_t i=0;i<in.size();++i) o[i]=in[i].vertex[2]; return o; }
inline Vec_d getMCParticleVertex_r(const Vec_MCPart& in){ Vec_d o(in.size()); for(size_t i=0;i<in.size();++i) o[i]=std::hypot(in[i].vertex[0],in[i].vertex[1]); return o; }

inline Vec_d getMCParticleVertex_theta(const Vec_MCPart& in,bool deg=false){
    Vec_d o(in.size());
    for(size_t i=0;i<in.size();++i){
        const auto& p=in[i];
        double r=std::sqrt(p.vertex[0]*p.vertex[0]+p.vertex[1]*p.vertex[1]+p.vertex[2]*p.vertex[2]);
        double t=r>0?std::acos(p.vertex[2]/r):0.; if(deg) t*=180./M_PI; o[i]=t;
    } return o;
}
inline Vec_d getMCParticleVertex_phi(const Vec_MCPart& in,bool deg=false){
    Vec_d o(in.size());
    for(size_t i=0;i<in.size();++i){ double p=std::atan2(in[i].vertex[1],in[i].vertex[0]); if(deg) p*=180./M_PI; o[i]=p; }
    return o;
}
inline Vec_d getMCParticleEnergy(const Vec_MCPart& in){
    Vec_d o(in.size());
    for(size_t i=0;i<in.size();++i){
        const auto& p=in[i];
        o[i]=std::sqrt(p.momentum[0]*p.momentum[0]+p.momentum[1]*p.momentum[1]+p.momentum[2]*p.momentum[2]+p.mass*p.mass);
    } return o;
}
inline Vec_b filterGenStatus1(const Vec_MCPart& in){ Vec_b o(in.size()); for(size_t i=0;i<in.size();++i) o[i]=(in[i].generatorStatus==1); return o; }
