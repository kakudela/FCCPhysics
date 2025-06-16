# Particle gun (HepMC3)
Generator of particle gun in the `hepmc3` format. Uniform in polar angle and exponential in momentum range.

## Compile
The gun needs to compiled first:

```
./install_gunHEPMC3.sh
```

## Generate events

Create an input file (`gun.input`):

    npart 1
    theta_range 90.0,90.0
    mom_range 20.0,80.0
    pid_list 13,-13
    nevents 100000


This generates 100000 events (`nevents`), each event containing 1 particle (`npart_range`) with PID chosen randomly from `pid_list `, at 90 degrees and within a momentum range of 20 to 80 GeV. To run, do:


    ./run_gunHEPMC3_singularity.sh gun.input


The output is stored in the `gun.hepmc` file.



# Use Delphes to run on HepMC3
The HepMC3 file can be used with Delphes to generate the detector response and store it in the EDM4HEP format. For this, we first need to install the HepMC reader for Delphes `DelphesHepMC_EDM4HEP`, which can be done using the installer script:

    chmod 755 install_DelphesHEPMC.sh
    ./install_DelphesHEPMC.sh

These commands have to be executed once. If you open a new terminal, you only need to load the proper environmental variables by doing:

    source env.sh

To run 

    DelphesHepMC_EDM4HEP ../delphes/cards/IDEA_2T.tcl ../delphes/cards/output.tcl mu_theta_10-90_p_5.root mu_theta_10-90_p_5.hepmc



