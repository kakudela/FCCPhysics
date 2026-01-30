import sys, os, glob, shutil
import time
import argparse
import logging
import subprocess
import random
import socket
import re

logging.basicConfig(format='%(levelname)s: %(message)s')
logger = logging.getLogger("fcclogger")
logger.setLevel(logging.INFO)

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--submit", action='store_true', help="Submit to batch system")
parser.add_argument("--gun_file", type=str, help="Input gun file")
parser.add_argument("--gp_dir", type=str, help="Input directory containing .pairs files from Guinea pig")
parser.add_argument("--geometry_file", type=str, help="Geometry file (locally)")
parser.add_argument("--geometry_tag", type=str, help="Geometry tag (from Key4hep)")
parser.add_argument("--condor_queue", type=str, help="Condor priority", choices=["espresso", "microcentury", "longlunch", "workday", "tomorrow", "testmatch", "nextweek"], default="longlunch")
parser.add_argument("--condor_priority", type=str, help="Condor priority", default="group_u_FCC.local_gen")
parser.add_argument("--storagedir", type=str, help="Base directory to save the samples", default="/ceph/submit/data/group/fcc/kudela/ddsim/guineapig_lattices/")  # --crossingAngleBoost 0.015
parser.add_argument("--use_storagedir", action="store_true", help="Store outputs under --storagedir (gp mode). Default: store outputs in --gp_dir.")
parser.add_argument("--cms_pool", action="store_true", help="Submit to CMS pool")
parser.add_argument("--osg_pool", action="store_true", help="Submit to OSG pool (Open Science Grid)")
parser.add_argument("--nJobs", help="Maximum number of events", type=int, default=9e99)
parser.add_argument("--maxMemory", help="Maximum job memory", type=float, default=2000)
parser.add_argument("-x", "--xrootd", action='store_true', help="Use XrootD transfer")
args = parser.parse_args()

# GP_STACK = "/cvmfs/sw.hsf.org/key4hep/setup.sh -r 2024-03-10" # CLD_o2_v05
GP_STACK = "/cvmfs/sw.hsf.org/key4hep/setup.sh -r 2025-05-29"  # CLD_o2_v07
# GP_STACK = "/cvmfs/sw.hsf.org/key4hep/setup.sh -r 2025-01-28" # IDEA_o1_v03
SINGULARITY = "/cvmfs/singularity.opensciencegrid.org/opensciencegrid/osgvo-el9:latest"
HOSTNAME = socket.gethostname()

# python submit.py --gp_dir /ceph/submit/data/group/fcc/kudela/guineapig_studies/FCCee_Z_GHC_V25p1_FCCee_Z128_0T_grids8_vtx000 --geometry_file CLD_o2_v07_2T/CLD_o2_v07_2T.xml --osg_pool

# python submit.py --gp_dir /ceph/submit/data/group/fcc/ee/detector/guineapig_studies/FCCee_Z_GHC_V25p1_FCCee_Z128_0T_grids8 --storagedir /ceph/submit/data/group/fcc/kudela/guineapig_studies/FCCee_Z_GHC_V25p1_FCCee_Z128_0T_grids8_2Tddsim/ --geometry_file CLD_o2_v07_2T/CLD_o2_v07_2T.xml --osg_pool


def get_voms_proxy_path():
    try:
        output = subprocess.check_output(['voms-proxy-info'], text=True)
        for line in output.splitlines():
            if line.strip().startswith('path'):
                return line.split(':', 1)[1].strip()
    except subprocess.CalledProcessError as e:
        print(f"Error running voms-proxy-info: {e}")
    return None

pdg_map = {
    13: 'mu-',
    -13: 'mu+',
    11: 'e-',
    -11: 'e+',
    22: 'gamma',
    211: 'pi+',
    -211: 'pi-',
}

class DDSimProducer:

    def __init__(self, storagedir, args):
        self.geometry_tag = None
        self.geometry_file = None
        if args.geometry_tag:
            self.geometry_tag = args.geometry_tag  # automatically expanded with env vars
            self.geometry_name = os.path.splitext(os.path.basename(self.geometry_tag))[0]
            if not os.path.exists(self.geometry_tag):
                logger.error(f"Geometry file {self.geometry_tag} does not exist, exiting")
                sys.exit(1)
        elif args.geometry_file:
            self.geometry_file = os.path.abspath(os.path.expandvars(args.geometry_file))
            self.geometry_dir = os.path.dirname(self.geometry_file)
            self.geometry_name = os.path.splitext(os.path.basename(self.geometry_file))[0]

            if not os.path.exists(self.geometry_file):
                logger.error(f"Geometry file {self.geometry_file} does not exist, exiting")
                sys.exit(1)

        self.storagedir = storagedir
        self.cwd = os.getcwd()

        logger.info(f"Found geometry {self.geometry_name}")

        self.local_dir = f"{self.cwd}/local_ddsim/{self.geometry_name}/"

        self.condor_queue = args.condor_queue
        self.condor_priority = args.condor_priority
        self.args = args
        self.ddsim_args = ""

    def generate_gp_submit(self):

        gun, gp = False, False
        seeds = []

        if args.gun_file and os.path.exists(args.gun_file):
            gun = True
            self.gun_file = os.path.abspath(args.gun_file)
            self.name = os.path.basename(os.path.normpath(self.gun_file)).split('.')[0]
            logger.info(f"Found gun file {self.name}")

            multiplicity, particle, nevents = None, None, None
            thetaMin, thetaMax, momentumMin, momentumMax = None, None, None, None
            with open(self.gun_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    key, value = line.split(None, 1)  # split at first whitespace
                    if key == 'npart':
                        multiplicity = value
                    if key == 'nevents':
                        nevents = value
                    if key == 'theta_range':
                        thetaMin, thetaMax = value.split(',')
                    if key == 'mom_range':
                        momentumMin, momentumMax = value.split(',')
                    if key == 'pid_list':
                        if ',' in value:
                            logger.error(f"ddsim does not support a list of particles {value}")
                            quit()
                        particle = pdg_map[int(value)]

            self.ddsim_args = (
                f"--enableGun --gun.particle {particle} --gun.multiplicity {multiplicity} "
                f"--gun.thetaMin '{thetaMin}*deg' --gun.thetaMax '{thetaMax}*deg' "
                f"--gun.momentumMin '{momentumMin}*GeV' --gun.momentumMax '{momentumMax}*GeV' "
                f"--random.seed ${{seed}} --numberOfEvents {nevents} --gun.distribution uniform"
            )

            # For gun mode, keep seeds as numeric random seeds
            seeds = [str(random.randint(10000, 99999)) for _ in range(int(args.nJobs))]

        elif args.gp_dir and os.path.exists(args.gp_dir):
            gp = True
            self.gp_dir = os.path.abspath(args.gp_dir)
            self.name = os.path.basename(os.path.normpath(self.gp_dir))
            logger.info(f"Found guineapig directory {self.name}")

            logger.info(f"Get .pairs files from {self.gp_dir}")
            input_files = glob.glob(f"{self.gp_dir}/*.pairs")
            if len(input_files) == 0:
                logger.error(f"Input directory {self.gp_dir} does not contain .pairs files, exiting")
                sys.exit(1)
            logger.info(f"Found {len(input_files)} pairs files")

            # Build job list supporting BOTH output_*.pairs and output0_*.pairs in same directory.
            # Each job row will be: (ID, INFILE, OUTFILE)
            jobs = []
            input_files = sorted(input_files)

            for n, f in enumerate(input_files):
                if n >= args.nJobs:
                    break

                infile = os.path.basename(f)  # e.g. "output_927450.pairs" or "output0_927450.pairs"
                base = infile.replace(".pairs", "")

                m = re.search(r"(\d+)$", base)
                if not m:
                    logger.error(f"Could not parse numeric ID from filename: {infile}")
                    sys.exit(1)
                idnum = m.group(1)

                if base.startswith("output0_"):
                    outfile = f"ddsim0_{idnum}.root"
                elif base.startswith("output_"):
                    outfile = f"ddsim_{idnum}.root"
                else:
                    # If you have any other .pairs in the folder, skip them safely.
                    logger.warning(f"Skipping unexpected .pairs name (not output_/output0_): {infile}")
                    continue

                jobs.append((idnum, infile, outfile))

            if len(jobs) == 0:
                logger.error("No valid output_*.pairs or output0_*.pairs found to submit.")
                sys.exit(1)

            # ddsim args for gp mode (input/output handled explicitly in job script via infile/outfile)
            self.ddsim_args = f"--numberOfEvents -1"

        else:
            logger.error("Please provide input guinea pig directory or a gun file")
            quit()

        # preparing submission dirs
        # Default: store outputs in gp_dir for gp mode
        # If --use_storagedir is given, store under storagedir
        if gp:
            if args.use_storagedir:
                self.out_dir = f"{self.storagedir}/{self.geometry_name}/{self.name}/"
                self.log_dir = f"{self.out_dir}/logs/"
            else:
                self.out_dir = self.gp_dir
                self.log_dir = f"{self.gp_dir}/logs/"
        else:
            # gun mode has no gp_dir, always use storagedir structure
            self.out_dir = f"{self.storagedir}/{self.geometry_name}/{self.name}/"
            self.log_dir = f"{self.out_dir}/logs/"

        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir)
        if not os.path.exists(self.log_dir):
            os.makedirs(self.log_dir)

        # pack the geometry files to be shipped
        if self.geometry_file is not None:
            logger.info("Packing geometry files")
            geometry_sandbox = f"{self.out_dir}/geometry.tar"
            os.system(f"cd {self.geometry_dir} && tar -cf {geometry_sandbox} .")
            compactFile = f"geometry/{self.geometry_name}.xml"
            geometry_action = "mkdir geometry && tar -xf geometry.tar -C geometry"
        else:
            compactFile = self.geometry_tag
            geometry_action = ""

        out_dir_xrd = self.out_dir.replace("/ceph/submit", "")

        if args.xrootd:
            xrootd_cfg = f"""
if [ -d /cvmfs/grid.cern.ch/etc/grid-security/certificates ]; then
    export X509_CERT_DIR=/cvmfs/grid.cern.ch/etc/grid-security/certificates
else
    echo "ERROR: Grid certs not available on CVMFS" >&2
    exit 1
fi

voms-proxy-info -all

xrdcp --version
xrdcp -d 3 "${{outfile}}" root://submit50.mit.edu/{out_dir_xrd}/"${{outfile}}"
rc=$?
if [ $rc -ne 0 ]; then
    echo "xrdcp failed with exit code $rc" >&2
fi
echo "Copy output file via XrootD done"
            """
        else:
            xrootd_cfg = ""

        # Job script:
        # - gun mode: seed only (numeric),
        # - gp mode: we pass (idnum, infile, outfile) and produce outfile
        job = f"""
#!/bin/bash
set -e

unset LD_LIBRARY_PATH
unset PYTHONHOME
unset PYTHONPATH

echo "Release:"
cat /proc/version

echo "Hostname:"
hostname

echo "List current working dir"
ls -lrt

# Args:
# gun: $1 = seed
# gp : $1 = idnum, $2 = infile, $3 = outfile
arg1="$1"
arg2="$2"
arg3="$3"

# unpack geometry
{geometry_action}

echo "Checking CVMFS..."
if ! timeout 15 ls /cvmfs/sw.hsf.org/key4hep/setup.sh >/dev/null 2>&1; then
    echo "ERROR: CVMFS or key4hep setup.sh not found!" >&2
    exit 1
fi

echo "CVMFS looks OK. Sourcing Key4Hep"
if ! source {GP_STACK} >/dev/null 2>&1; then
    echo "ERROR: Failed to source Key4Hep setup" >&2
    exit 1
fi
echo "Key4Hep sourcing successful"

echo "Start ddsim"
SECONDS=0

if [ -z "$arg2" ] || [ -z "$arg3" ]; then
    # gun mode
    export seed="$arg1"
    echo "Gun seed: $seed"
    ddsim --compactFile {compactFile} --outputFile "${{seed}}_sim.root" {self.ddsim_args}
    outfile="${{seed}}_sim.root"
else
    # gp mode
    idnum="$arg1"
    infile="$arg2"
    outfile="$arg3"
    echo "GP id: $idnum"
    echo "GP infile: $infile"
    echo "GP outfile: $outfile"
    ddsim --compactFile {compactFile} --outputFile "${{outfile}}" --inputFiles "${{infile}}" {self.ddsim_args}
fi

duration=$SECONDS

rc=$?
if [ $rc -ne 0 ]; then
    echo "ddsim failed with exit code $rc" >&2
    exit $rc
fi

if [ ! -f "${{outfile}}" ]; then
    echo "ddsim did not produce ${{outfile}}" >&2
    exit 1
fi

echo "Done ddsim"
ls -lrt

if ! rootls -t "${{outfile}}" >/dev/null 2>&1; then
  echo "File got corrupted" >&2
  exit 1
fi

{xrootd_cfg}

echo "Duration: $(($duration))
"
        """

        # make executable script
        submitFn = f"{self.out_dir}/run_ddsim.sh"
        with open(submitFn, "w") as fOut:
            fOut.write(job)
        subprocess.getstatusoutput(f"chmod 777 {submitFn}")

        # make condor submission script
        condorFn = f"{self.out_dir}/condor.cfg"
        fOut = open(condorFn, 'w')

        fOut.write('universe       = vanilla\n')
        fOut.write(f'initialdir     = {self.out_dir}\n')
        fOut.write(f'executable     = {submitFn}\n')

        # Arguments differ between gp and gun
        if gp:
            fOut.write('arguments      = $(IDNUM) $(INFILE) $(OUTFILE)\n')
        else:
            fOut.write('arguments      = $(SEED)\n')

        fOut.write('Log            = logs/condor_job.$(ClusterId).$(ProcId).log\n')
        fOut.write('Output         = logs/condor_job.$(ClusterId).$(ProcId).out\n')
        fOut.write('Error          = logs/condor_job.$(ClusterId).$(ProcId).error\n')

        fOut.write('should_transfer_files = YES\n')
        fOut.write('when_to_transfer_output = ON_EXIT\n')

        # Inputs
        if gp:
            if self.geometry_file is not None:
                # Transfer the exact .pairs filename (works for output_ and output0_ in same directory)
                fOut.write(f'transfer_input_files = {self.gp_dir}/$(INFILE),{geometry_sandbox}\n')
            else:
                fOut.write(f'transfer_input_files = {self.gp_dir}/$(INFILE)\n')
        else:
            if self.geometry_file is not None:
                fOut.write(f'transfer_input_files = {geometry_sandbox}\n')

        # Outputs
        if args.xrootd:
            fOut.write('transfer_output_files = ""\n')
        else:
            if gp:
                fOut.write('transfer_output_files = $(OUTFILE)\n')
            else:
                fOut.write('transfer_output_files = $(SEED)_sim.root\n')

        fOut.write('on_exit_remove = (ExitBySignal == False) && (ExitCode == 0)\n')
        fOut.write('max_retries    = 150\n')
        fOut.write('periodic_release = (HoldReasonCode == 12 && NumJobStarts < 150)\n')
        fOut.write('on_exit_hold = False\n')
        fOut.write(f'RequestMemory  = {args.maxMemory}\n')

        if args.xrootd:
            proxy_path = get_voms_proxy_path()
            if proxy_path:
                os.system(f"cp {proxy_path} {self.log_dir}/")
                fOut.write('use_x509userproxy     = True\n')
                fOut.write(f'x509userproxy         = logs/{os.path.basename(proxy_path)}\n')

        # global pool
        if args.cms_pool:
            fOut.write('+DESIRED_Sites = "T2_AT_Vienna,T2_BE_IIHE,T2_BE_UCL,T2_BR_SPRACE,T2_BR_UERJ,T2_CH_CERN,T2_CH_CERN_AI,T2_CH_CERN_HLT,T2_CH_CERN_Wigner,T2_CH_CSCS,T2_CH_CSCS_HPC,T2_CN_Beijing,T2_DE_DESY,T2_DE_RWTH,T2_EE_Estonia,T2_ES_CIEMAT,T2_ES_IFCA,T2_FI_HIP,T2_FR_CCIN2P3,T2_FR_GRIF_IRFU,T2_FR_GRIF_LLR,T2_FR_IPHC,T2_GR_Ioannina,T2_HU_Budapest,T2_IN_TIFR,T2_IT_Bari,T2_IT_Legnaro,T2_IT_Pisa,T2_IT_Rome,T2_KR_KISTI,T2_MY_SIFIR,T2_MY_UPM_BIRUNI,T2_PK_NCP,T2_PL_Swierk,T2_PL_Warsaw,T2_PT_NCG_Lisbon,T2_RU_IHEP,T2_RU_INR,T2_RU_ITEP,T2_RU_JINR,T2_RU_PNPI,T2_RU_SINP,T2_TH_CUNSTDA,T2_TR_METU,T2_TW_NCHC,T2_UA_KIPT,T2_UK_London_IC,T2_UK_SGrid_Bristol,T2_UK_SGrid_RALPP,T2_US_Caltech,T2_US_Florida,T2_US_MIT,T2_US_Nebraska,T2_US_Purdue,T2_US_UCSD,T2_US_Vanderbilt,T2_US_Wisconsin,T3_CH_CERN_CAF,T3_CH_CERN_DOMA,T3_CH_CERN_HelixNebula,T3_CH_CERN_HelixNebula_REHA,T3_CH_CMSAtHome,T3_CH_Volunteer,T3_US_HEPCloud,T3_US_NERSC,T3_US_OSG,T3_US_PSC,T3_US_SDSC"\n')
            fOut.write('+AccountingGroup      = "analysis.jaeyserm"\n')
            fOut.write('+SingularityImage       = "/cvmfs/singularity.opensciencegrid.org/opensciencegrid/osgvo-el9:latest"\n')
            fOut.write('+ProjectName            = "MIT_submit"\n')
            fOut.write('+SingularityBindCVMFS   = True\n')
            fOut.write('Requirements          = ( BOSCOCluster =!= "t3serv008.mit.edu" && BOSCOCluster =!= "ce03.cmsaf.mit.edu" && BOSCOCluster =!= "eofe8.mit.edu")\n')

        # OSG pool
        elif args.osg_pool:
            fOut.write(f'+SingularityImage       = "{SINGULARITY}"\n')
            fOut.write('+SINGULARITY_BIND_EXPR       = "/cvmfs"\n')
            fOut.write('+ProjectName            = "MIT_submit"\n')
            fOut.write('+SingularityBindCVMFS   = True\n')
            fOut.write(f'Requirements          = ( OSGVO_OS_STRING == "RHEL 9" && HAS_CVMFS_singularity_opensciencegrid_org == TRUE && HAS_SINGULARITY == TRUE &&  (GLIDEIN_Site == "Wisconsin" || GLIDEIN_Site == "UChicago")  )\n')

        elif 'mit.edu' in HOSTNAME:
            proxy_path = get_voms_proxy_path()
            if proxy_path:
                os.system(f"cp {proxy_path} {self.log_dir}/x509")
                fOut.write('use_x509userproxy     = True\n')
                fOut.write('x509userproxy         = logs/x509\n')
            fOut.write('Requirements          = ( BOSCOCluster =!= "t3serv008.mit.edu" && BOSCOCluster =!= "ce03.cmsaf.mit.edu" && BOSCOCluster =!= "eofe8.mit.edu")\n')
            fOut.write('+DESIRED_Sites = "mit_tier2,mit_tier3"\n')
            fOut.write(f'+SingularityImage       = "{SINGULARITY}"\n')
            fOut.write('+SINGULARITY_BIND_EXPR       = "/cvmfs,/etc/grid-security"\n')
            fOut.write('+SingularityBindCVMFS   = True\n')
            fOut.write('Requirements          = ( BOSCOCluster =!= "t3serv008.mit.edu" && BOSCOCluster =!= "ce03.cmsaf.mit.edu" && BOSCOCluster =!= "eofe8.mit.edu")\n')

        elif 'cern.ch' in HOSTNAME:
            fOut.write(f'+JobFlavour    = "{self.condor_queue}"\n')
            fOut.write(f'+AccountingGroup = "{self.condor_priority}"\n')

        # Queue
        if gp:
            queue_txt = os.path.join(self.out_dir, "queue.txt")
            with open(queue_txt, "w") as qf:
                for (idnum, infile, outfile) in jobs:
                    qf.write(f"{idnum} {infile} {outfile}\n")
            fOut.write(f'queue IDNUM,INFILE,OUTFILE from {queue_txt}\n')
        else:
            seedsStr = ' \n '.join([str(s) for s in seeds])
            fOut.write(f'queue SEED in ( \n {seedsStr} \n)\n')

        fOut.close()

        subprocess.getstatusoutput(f"chmod 777 {condorFn}")
        os.system(f"condor_submit {condorFn}")

        logger.info(f"Written to {self.out_dir}")

def main():
    producer = DDSimProducer(args.storagedir, args)
    producer.generate_gp_submit()

if __name__ == "__main__":
    main()
