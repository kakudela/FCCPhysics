

import sys, os, glob
import argparse
import logging
import subprocess

logging.basicConfig(format='%(levelname)s: %(message)s')
logger = logging.getLogger("fcclogger")
logger.setLevel(logging.INFO)


parser = argparse.ArgumentParser()
parser.add_argument("--gen_card", type=str, help="Generator card", default="cards/Zmumu_ecm240.cmd")
parser.add_argument("--delphes_card", type=str, help="Detector card", default="cards/IDEA_2T.tcl")
parser.add_argument("--steering_card", type=str, help="Steering card", default="cards/output.tcl")
parser.add_argument("--output_dir", type=str, help="Step", default="samples")


args = parser.parse_args()

def main():

    STACK = "/cvmfs/sw.hsf.org/spackages7/key4hep-stack/2023-04-08/x86_64-centos7-gcc11.2.0-opt/urwcv/setup.sh"
    #STACK = "/cvmfs/sw.hsf.org/spackages6/key4hep-stack/2022-12-23/x86_64-centos7-gcc11.2.0-opt/ll3gi/setup.sh"
    STACK = "/cvmfs/sw.hsf.org/key4hep/setup.sh -r 2023-11-23"

    gen_card = args.gen_card
    delphes_card = args.delphes_card
    steering_card = args.steering_card
    output_dir = args.output_dir

    if not os.path.exists(gen_card):
        logger.error(f"Generator card {gen_card} does not exist.")
        quit()
    if not os.path.exists(delphes_card):
        logger.error(f"Delphes card {delphes_card} does not exist.")
        quit()
    if not os.path.exists(steering_card):
        logger.error(f"Steering card {steering_card} does not exist.")
        quit()
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    detector = os.path.splitext(os.path.basename(delphes_card))[0]
    process = os.path.splitext(os.path.basename(gen_card))[0]
    logger.info(f"Running process {process} with detector {detector}")
    output_file = f"{output_dir}/{detector}_{process}.root"
    if os.path.exists(output_file):
        logger.error(f"Output file {output_file} already exists, please remove.")
        quit()

    logger.info(f"Start generating events")
    cmd = f"source {STACK} && DelphesPythia8_EDM4HEP {delphes_card} {steering_card} {gen_card} {output_file}"
    subprocess.run(cmd, shell=True, env={})
    logger.info(f"Done! Output saved to {output_file}")

if __name__ == "__main__":
    main()
