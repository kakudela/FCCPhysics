set -euo pipefail

# where to store all outputs and the timing log
outbase="/ceph/submit/data/group/fcc/kudela/gp_timing_studies"
log="${outbase}/timings.csv"

# guinea binary
guinea="/ceph/submit/data/group/fcc/ee/guineapig/guinea"
# blocks in acc.dat
acc_block="FCCee_Z_GHC_V25p1"
par_block="TEST_PAR"

mkdir -p "${outbase}"

# write CSV header once if file doesn't exist
if [ ! -f "${log}" ]; then
  echo "timestamp,host,grids,track_pairs,output_dir,real_s,user_s,sys_s,maxrss_kb" > "${log}"
fi

# function to update grids and track_pairs inside acc.dat (in-place)
set_param () {
  local key="$1"
  local val="$2"
  # replaces lines like "grids=8 ;" or "track_pairs=2;"
  sed -i -E "s/^([[:space:]]*${key}[[:space:]]*=[[:space:]]*).*/\1${val};/" acc.dat
}

run_one () {
  local grids="$1"
  local tp="$2"
  local name="$3"

  set_param "grids" "${grids}"
  set_param "track_pairs" "${tp}"

  # make sure seed does not advance between runs
  rm -f rndm.save

  local outdir="${outbase}/${name}"
  mkdir -p "${outdir}"

  local ts
  ts="$(date +%Y-%m-%dT%H:%M:%S)"
  local host
  host="$(hostname)"

  # run and capture timing to a temp file
  local tfile
  tfile="$(mktemp)"
  /usr/bin/time -f "%e,%U,%S,%M" -o "${tfile}" \
  bash -c "
    cd '${outdir}' && \
    '${guinea}' --acc_file '${PWD}/acc.dat' '${acc_block}' '${par_block}' .
  "


  # parse timing fields
  local real user sys rss
  IFS=',' read -r real user sys rss < "${tfile}"
  rm -f "${tfile}"

  echo "${ts},${host},${grids},${tp},${outdir},${real},${user},${sys},${rss}" >> "${LOG}"
  echo "logged: grids=${grids} track_pairs=${tp} out=${outdir} real_s=${real}"
}

# scan: grids number, track pairs, output subdir name
run_one 1 0 "FCCee_Z_GHC_V25p1_2T_grids1_noTracking"
run_one 1 1 "FCCee_Z_GHC_V25p1_2T_grids1_tracking"
run_one 2 0 "FCCee_Z_GHC_V25p1_2T_grids2_noTracking"
run_one 2 1 "FCCee_Z_GHC_V25p1_2T_grids2_tracking"
run_one 3 0 "FCCee_Z_GHC_V25p1_2T_grids3_noTracking"
run_one 3 1 "FCCee_Z_GHC_V25p1_2T_grids3_tracking"
run_one 4 0 "FCCee_Z_GHC_V25p1_2T_grids4_noTracking"
run_one 4 1 "FCCee_Z_GHC_V25p1_2T_grids4_tracking"
run_one 5 0 "FCCee_Z_GHC_V25p1_2T_grids5_noTracking"
run_one 5 1 "FCCee_Z_GHC_V25p1_2T_grids5_tracking"
run_one 6 0 "FCCee_Z_GHC_V25p1_2T_grids6_noTracking"
run_one 6 1 "FCCee_Z_GHC_V25p1_2T_grids6_tracking"
run_one 7 0 "FCCee_Z_GHC_V25p1_2T_grids7_noTracking"
run_one 7 1 "FCCee_Z_GHC_V25p1_2T_grids7_tracking"
run_one 8 0 "FCCee_Z_GHC_V25p1_2T_grids8_noTracking"
run_one 8 1 "FCCee_Z_GHC_V25p1_2T_grids8_tracking"
