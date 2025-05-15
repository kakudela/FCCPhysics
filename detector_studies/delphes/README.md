# Detector Studies using DELPHES
This repository contains a basic setup to generate local Delphes events given a physics process with Pythia, and analyze the output to extract the detector performance.

Environment:

	source /work/submit/jaeyserm/software/FCCAnalyses/setup.sh


#### Generate events

To generate events we need:

- A physics process (e.g. Z decay to 2 muons)
- A detector description (e.g. IDEA with gaseous tracker or CLD with silicon tracker).

Both the physiscs process and the detector description are defined in input cards (text files), in the `cards` directory:

- `Zmumu_ecm240.cmd` Z decay to 2 muons at center-of-mass energy of 240 GeV (Pythia card)
- `IDEA_2T.tcl` IDEA detector with 2T magnetic field (Delphes card)
- `CLD_2T.tcl` CLD detector with 2T magnetic field (Delphes card)

To generate events, do:

	python generate_local.py --delphes_card cards/CLD_2T.tcl --gen_card cards/Zmumu_ecm240.cmd

It generates N events where N is defined in the Pythia card. It stores a single file with all events in the `output` folder.

#### Analyze events

We extract the momentum resolution, defined as `(p_reco p_gen)/p_gen` using the following script:

	python analysis_resolution.py

This script runs over all the events generated above (specify the process/detector in the `analysis_resolution.py`), and creates another file containing the histogram of the resolution. To plot it and extract the resolution, execute:

	python plot_resolution.py

The output of this script creates the resolution plot, where the resolution is displayed.