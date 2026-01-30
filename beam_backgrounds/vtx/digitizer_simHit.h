#pragma once
#include "functions.h"
#include <math.h>

// Jan's digitizer
// detailed explanation of this script in README.md, please check it out!

struct PixelDigi {
    int row;
    int col;
    float adc; // MeV stored here
};

struct BarrelGeometry {
    float R;  // layer radius, mm
    float z;  // half-length, mm
    float pX = 25; // pitch along rphi, um
    float pY = 25; // pitch along z,   um
    float c;  // circumference, mm
    BarrelGeometry() {}
    BarrelGeometry(float R_, float z_, float pX_, float pY_) : R(R_), z(z_), pX(pX_), pY(pY_), c(2*M_PI*R_) {}
};

struct EndcapGeometry {
    float rMax;  // half-size of square canvas containing the disk, mm
    float zLo;   // inclusive lower bound on |z|, mm
    float zHi;   // exclusive upper bound on |z|, mm
    float pX = 25; // pitch along x, um
    float pY = 25; // pitch along y, um
    float w;    // canvas width,  mm
    float h;    // canvas height, mm
    EndcapGeometry() {}
    EndcapGeometry(float rMax_, float zLo_, float zHi_, float pX_, float pY_)
        : rMax(rMax_), zLo(zLo_), zHi(zHi_), pX(pX_), pY(pY_), w(2*rMax_), h(2*rMax_) {}
};

using Vec_PixelDigi = ROOT::VecOps::RVec<PixelDigi>;

// Barrel: unroll to (r*phi, z)
Vec_PixelDigi DigitizerSimHitBarrel(const BarrelGeometry& geo, Vec_f hits_x, Vec_f hits_y, Vec_f hits_z, Vec_f hits_e, float thrs) {

    Vec_PixelDigi digis;

    for(int i=0; i<(int)hits_e.size(); i++) {
        float depositEnergy = hits_e[i] * 1000; // MeV
        if(depositEnergy < thrs) continue;

        float phi = std::atan2(hits_y[i], hits_x[i]);
        phi = (phi > 0 ? phi : (2*M_PI + phi)); // 0..2pi
        float x_pos = geo.R*phi*1000.f; // um
        float y_pos = hits_z[i]*1000.f; // um

        int row = static_cast<int>(x_pos / geo.pX);
        int col = static_cast<int>(y_pos / geo.pY);

        digis.push_back({row, col, depositEnergy});
    }

    return digis;
}

// End-cap: planar x-y projection, both sides selected by |z| (wip)
Vec_PixelDigi DigitizerSimHitEndcap(const EndcapGeometry& geo, Vec_f hits_x, Vec_f hits_y, Vec_f hits_z, Vec_f hits_e, float thrs) {

    Vec_PixelDigi digis;

    const float x_off_um = geo.rMax * 1000.f; // [-rMax,+rMax] -> [0,2*rMax]
    const float y_off_um = geo.rMax * 1000.f;

    for(int i=0; i<(int)hits_e.size(); i++) {
        const float az = std::fabs(hits_z[i]); // use both +/- z
        if(az < geo.zLo || az >= geo.zHi) continue;

        float depositEnergy = hits_e[i] * 1000.f; // MeV
        if(depositEnergy < thrs) continue;

        float u_um = hits_x[i]*1000.f + x_off_um;
        float v_um = hits_y[i]*1000.f + y_off_um;

        int row = static_cast<int>(u_um / geo.pX);
        int col = static_cast<int>(v_um / geo.pY);

        digis.push_back({row, col, depositEnergy});
    }

    return digis;
}

// Barrel occupancy on unrolled rectangle
Vec_f getOccupancyBarrel(const BarrelGeometry& geo, Vec_PixelDigi digis, int nR, int nC, bool npx=false) {

    int n_row_pixels_per_window, n_col_pixels_per_window;
    if(npx) {
        n_row_pixels_per_window = nR;
        n_col_pixels_per_window = nC;
    }
    else {
        n_row_pixels_per_window = std::ceil(geo.c * 1000.f / nR / geo.pX);
        n_col_pixels_per_window = std::ceil(geo.z * 1000.f * 2.0f / nC / geo.pY);
    }
    int n_col_pixels_half_z = geo.z * 1000.f / geo.pY;

    Vec_i hitmap = Vec_i(nR*nC, 0);
    for(const auto& b : digis) {
        int r = std::floor(b.row/n_row_pixels_per_window);
        int c = std::floor((b.col + n_col_pixels_half_z)/n_col_pixels_per_window);
        if(r>=0 && r<nR && c>=0 && c<nC) hitmap[r*nC+c] += 1;
    }

    int maxCount = hitmap.empty() ? 0 : *std::max_element(hitmap.begin(), hitmap.end());
    float sum = hitmap.empty() ? 0.f : std::accumulate(hitmap.begin(), hitmap.end(), 0.0f);
    float avgCount = (hitmap.empty() ? 0.f : sum / (float)hitmap.size());
    int totCount = (int)digis.size();

    Vec_f ret;
    ret.push_back(maxCount);
    ret.push_back(avgCount);
    ret.push_back(totCount);
    return ret;
}

// Planar occupancy on x-y canvas of size (2*rMax) x (2*rMax)
Vec_f getOccupancyPlanar(const EndcapGeometry& geo, Vec_PixelDigi digis, int nR, int nC, bool npx=false) {

    int n_row_pixels_per_window, n_col_pixels_per_window;
    if(npx) {
        n_row_pixels_per_window = nR;
        n_col_pixels_per_window = nC;
    }
    else {
        n_row_pixels_per_window = std::ceil(geo.w * 1000.f / nR / geo.pX);
        n_col_pixels_per_window = std::ceil(geo.h * 1000.f / nC / geo.pY);
    }

    Vec_i hitmap = Vec_i(nR*nC, 0);
    for(const auto& b : digis) {
        int r = std::floor(b.row/n_row_pixels_per_window);
        int c = std::floor(b.col/n_col_pixels_per_window);
        if(r>=0 && r<nR && c>=0 && c<nC) hitmap[r*nC+c] += 1;
    }

    int maxCount = hitmap.empty() ? 0 : *std::max_element(hitmap.begin(), hitmap.end());
    float sum = hitmap.empty() ? 0.f : std::accumulate(hitmap.begin(), hitmap.end(), 0.0f);
    float avgCount = (hitmap.empty() ? 0.f : sum / (float)hitmap.size());
    int totCount = (int)digis.size();

    Vec_f ret;
    ret.push_back(maxCount);
    ret.push_back(avgCount);
    ret.push_back(totCount);
    return ret;
}
