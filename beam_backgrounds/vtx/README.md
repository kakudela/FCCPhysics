___________________________________________________________________________________
___________________________________________________________________________________

DDSIM ANALYSIS STUFF!

Detector: CLD_o2_v07.
source "/cvmfs/sw.hsf.org/key4hep/setup.sh -r 2025-05-29"

___________________________________________________________________________________
___________________________________________________________________________________


analysis.py

Runs analysis for DDSim outputs, producing VTX + GEN histograms and occupancy numbers.

- Input expectation:
  The input directory contains DDSim ROOT files:
    * ddsim_*.root = tracked GuineaPig “pairs” stream
    * ddsim0_*.root = untracked GuineaPig “pairs0” stream
  analysis.py runs a separate full analysis for each stream that exists (pairs and/or pairs0).

- Core workflow:
  -> Discover + sanity-check input files
     - Globs ddsim_*.root and ddsim0_*.root
     - Uses functions.py file pruning helpers to skip zombies / missing trees
  -> Build RDataFrame over the "events" tree
     - Produces derived columns (hit kinematics, MC kinematics, selections)
     - Fills TH1/TH2 histograms (UNNORMALIZED; i.e. raw counts across all events)
  -> Digitize L1 barrel simhits (geometric pixel mapping)
     - Calls digitizer_simHit.h routines to map SimTrackerHits -> PixelDigi(row,col,adc)
     - Counts digis per event and per local window (for occupancy proxies)
  -> Compute scalar summary numbers
     - n_events and per-event averages (avg simhits, avg unique MC, etc.)
     - global occupancy and local occupancy (avg and p95; written as scalars)
     - fractions such as “frac_mc_sel_with_l1hit”
  -> Persist outputs
     - Writes one ROOT output file per stream (pairs / pairs0) containing:
       * all TH1/TH2 histograms (raw/unscaled)
       * scalar parameters as TParameter<...> (and small text summaries)
     - Writes PNG plots into the corresponding public_html directory structure and
       auto-copies index.php into each directory via ensure_public_dir(...) / ensure_index_php(...)

- Plotting behavior:
  - Default mode: fills + writes ROOT outputs AND produces PNGs
  - --plot: plot-only mode (re-reads the already-produced ROOT file and just makes PNGs)

- Magnetic field handling (GEN theta–pT acceptance plot):
  - The script contains a mapping from dataset/directory tag -> assumed B field for
    GEN acceptance overlays.
  - If run on a single directory without a known tag, it assumes B = 2T for the
    acceptance box definition.

- Output locations:
  - ROOT output file is written into the same directory tree as the PNGs (per stream).
  - VTX plots go under:
      /home/submit/kudela/public_html/fccee/beam_background/ddsim/vtx/<dataset>/<stream>/
    GEN plots go under:
      /home/submit/kudela/public_html/fccee/beam_background/ddsim/gen/<dataset>/<stream>/
  - The canonical ROOT output name pattern is:
      <dataset>_ddsim_analysis.root

- Key point for downstream scripts:
  - Histograms are NOT per-event normalized in analysis.py.
    compare.py and other consumers must divide by n_events when they want per-event rates.
  - Scalar summaries are stored as TParameter and are read by occupancy_studies.py.


___________________________________________________________________________________
___________________________________________________________________________________


compare.py

Overlays and compares 1D histograms produced by analysis.py across multiple samples.

- Inputs:
  - A group = { samples[], plots[], outdir }
  - Each sample points to an analysis.py ROOT output file

- Key behavior:
  - analysis.py histograms are raw (not per-event). compare.py always divides by n_events
    (read from TParameter<int>("n_events") in the ROOT file; fallback to 1).
  - Two normalization modes:
      * per_event: only divide by n_events
      * unit_area: divide by n_events AND normalize histogram to unit integral

- Plotting:
  - Clones hists out of the file (TH1.AddDirectory(False) behavior)
  - Applies consistent styling + optional legend placement overrides
  - Optional x_range trimming for view and for ymax estimation
  - Saves PNGs to:
      /home/submit/kudela/public_html/fccee/beam_background/ddsim/compare
    and ensures index.php is present.


___________________________________________________________________________________
___________________________________________________________________________________


occupancy_studies.py

Reads scalar parameters written by analysis.py and makes “parameter vs x” plots for
multiple independent studies (B-field scans, grid scans, granularity scans, etc).

- Inputs:
  - A study is defined by a list of points: (x value, ROOT file path)
  - ROOT files are analysis.py outputs (…_ddsim_analysis.root)

- Expected ROOT contents (scalars):
  Scalars are stored as TParameter<...> with exact names like:
    avg_simhits_barrel_l1
    avg_unique_mc_barrel_l1
    avg_simhits_endcap_l1
    avg_unique_mc_endcap_l1
    avg_global_occupancy_u6
    p95_global_occupancy_u6
    frac_mc_sel_with_l1hit

- Reader behavior:
  - Reads TParameter via GetVal()
  - Fallback: supports TNamed with a numeric Title/Name
  - If a point is missing a parameter, that point is skipped for that plot.


___________________________________________________________________________________
___________________________________________________________________________________


functions.py

Python-side utilities used across analysis.py, occupancy_studies.py, and compare.py.

- Directory / web output helpers:
  - ensure_dir: mkdir -p
  - ensure_index_php: copies index.php into a directory once
  - ensure_public_dir: constructs nested public_html paths and ensures index.php exists

- Histogram saving helpers:
  - save_hist1d: clones hist, optional scaling, axis labels, optional log-y, x-buffering
  - save_hist2d: clones hist, optional scaling, axis labels, buffers x/y, draws COLZ,
    forces no stats box for 2D plots
  - save_hist2d_with_box: same as save_hist2d but overlays a red TBox in user coords

- Occupancy math helpers:
  - sensor_area_mm2, pixel_area_mm2
  - hist_mean_and_p95
  - compute_occ_from_hist: avg/p95 hits and converts → occupancy proxy
  - compute_global_occ_from_hits: hit count → occupancy proxy

- Text output helper:
  - write_text: writes per-run summary text files

- File health / pruning:
  - filter_good_files: opens each ROOT file and checks tree existence, drops zombies
  - weed_files_on_open: attempts to build an RDataFrame; if it fails, removes the
    offending file (or trims) to avoid analysis crashes on truncated files


___________________________________________________________________________________
___________________________________________________________________________________


functions.h

Serves as a lightweight C++ helper library for ROOT RDataFrame, providing compiled
accessors and small algorithms used heavily by analysis.py.

- Defines RVec-based type aliases (Vec_f, Vec_i, Vec_b, etc.) and edm4hep
  collection aliases for SimTrackerHitData and MCParticleData, allowing concise
  RDataFrame Define(...) expressions.

- Provides fast getters for SimTrackerHit quantities:
  hit positions (x,y,z,r,phi,theta), energy deposition, momentum components, and
  secondary/primary tagging via the simhit quality bits.

- Provides MCParticle helpers:
  vertex coordinates, momentum components, pt, total momentum, energy, beta,
  generator status, and PDG ID (optionally abs-valued). These are used for MC
  selection, kinematics, and truth-level plots.

- Contains small vector utilities used inside RDataFrame: absolute value, cos(theta),
  log10 transforms, angle wrapping to [0,360),
  and safe elementwise division.

- Implements MC–hit bookkeeping tools:
  counting unique MCs with ≥1 hit, per-MC hit multiplicity, and consecutive-hit
  delta calculations (dz, dphi) for studying hit correlations.

- Implements consecutive-hit merging (the “dedup25um” logic):
  hits from the same MC particle that lie within a configurable spatial radius
  (25 um in analysis.py) are merged to suppress artificial hit splitting.

- Includes geometry-derived helpers such as the barrel incidence angle, computed
  relative to the local radial normal.

All functions are header-only and are JIT-compiled once via
ROOT.gInterpreter.Declare() at runtime.


___________________________________________________________________________________
___________________________________________________________________________________


digitizer_simHit.h

Implements a minimal, geometry-based digitization layer that converts selected
SimTrackerHits into pixel-level “digis” and derives hit-density information used
to estimate detector occupancy in analysis.py.

IMPORTANT: this is NOT a full detector response simulation.
It is a simplified geometric mapping intended for fast occupancy studies.

-----------------------------------------------
What a “digi” means here

A digi is represented by:

    struct PixelDigi {
        int row;
        int col;
        float adc;   // stored as deposited energy in MeV
    };

- row, col:
  integer pixel indices on an abstract pixel grid
  (unrolled barrel surface or planar endcap canvas).

- adc:
  simply the SimTrackerHit deposited energy converted to MeV.
  There is no electronics response, noise, shaping, or ADC modeling.

-----------------------------------------------
Geometry parameterization

Two simplified sensor geometries are defined:

BarrelGeometry:
- R (mm): barrel layer radius (L1 inner barrel uses 14.0 mm)
- z (mm): half-length of the barrel (109.5 mm for L1 inner)
- pX, pY (um): pixel pitch along r–phi and z directions
- c = 2*pi*R: circumference, used for unrolling

EndcapGeometry:
- rMax (mm): half-size of square canvas bounding the endcap disk
- zLo, zHi (mm): |z| selection window
- pX, pY (um): pixel pitch in x and y
- w, h = 2*rMax: canvas dimensions

These geometries define only how hits are mapped to pixel indices;
they do not encode material, thickness, or response effects.

-----------------------------------------------
SimHit -> PixelDigi mapping

Barrel digitization:
    DigitizerSimHitBarrel(...)

Inputs (from analysis.py):
- x, y, z hit positions (mm), already filtered to L1 inner barrel
- hit energy deposition (GeV)
- pixel pitch and barrel geometry
- an energy threshold (currently set to 0)

Procedure for each SimTrackerHit:
1. Convert deposited energy:
       E_mev = eDep_GeV * 1000
2. Apply threshold cut:
       if E_mev < threshold: skip
3. Compute azimuthal angle:
       phi = atan2(y, x) mapped to [0, 2*pi)
4. Unroll cylindrical surface:
       x_unrolled = R * phi  (mm) → converted to um
       y_unrolled = z (mm) → converted to um
5. Compute pixel indices:
       row = floor(x_unrolled / pX)
       col = floor(y_unrolled / pY)
6. Store PixelDigi(row, col, E_mev)

Each accepted SimTrackerHit produces exactly one digi.

Endcap digitization:
    DigitizerSimHitEndcap(...)

Differences from barrel:
- Hit must satisfy zLo ≤ |z| < zHi
- Hits are mapped directly onto a 2D x–y plane
- A fixed offset maps [-rMax, +rMax] mm -> [0, 2*rMax] mm before pixelization
- row/col are computed from x/y instead of r–phi/z

-----------------------------------------------
Connection to analysis.py (global hits)

In analysis.py, digitization is applied only to L1 inner barrel hits:

    b_digis = DigitizerSimHitBarrel(...)
    global_tot_hits = b_digis.size()

Thus:
- global_tot_hits = number of digis per event
- this is the raw input for the “global occupancy” calculation

analysis.py then converts this count into an occupancy proxy by:
- normalizing by sensor area
- scaling by pixel area
- inflating by cluster_size and safety_factor

The digitizer itself does not compute occupancies — it provides the
hit counts on which analysis.py operates.

-----------------------------------------------
Local occupancy via windowed hit maps

Local hit density is computed using:

    getOccupancyBarrel(...)
    getOccupancyPlanar(...)

These functions:
- divide the sensor into a grid of windows
- count how many digis fall into each window

Inputs:
- Vec<PixelDigi> digis
- nR, nC: either window size in pixels or number of windows
- npx flag:
    * true  -> nR x nC pixels per window
    * false -> nR x nC windows across sensor

For each window, a hit count is accumulated.

Return value (Vec_f of size 3):
- [0] = maximum hit count in any window (local hotspot)
- [1] = average hits per window
- [2] = total number of digis

-----------------------------------------------
How local occupancy is formed in analysis.py

analysis.py uses only the *maximum* window hit count:

    local_max_hits = occ_cfg[0]

This is converted to an occupancy estimate via:

- Define the physical window area
  (pixel-based or geometry-based)
- Convert hits -> hits/mm²
- Multiply by:
    * pixel area
    * cluster_size (effective pixels per hit)
    * safety_factor (engineering margin)
- Scale by 1e6 for “u6” plots

This same scaling model is used for both global and local occupancy so
they are directly comparable.

-----------------------------------------------
Key simplifying assumptions

- One SimTrackerHit -> one digi (no charge sharing)
- No time structure, integration window, or pileup model
- No electronic thresholding beyond a simple energy cut
- No clustering algorithm (cluster_size applied later as a scalar)
- No diffusion, Lorentz drift, or sensor effects

This digitizer is intentionally minimal and fast, designed to study
relative occupancy trends as a function of beam background conditions,
granularity, and magnetic field — not detailed detector response.


___________________________________________________________________________________
___________________________________________________________________________________
