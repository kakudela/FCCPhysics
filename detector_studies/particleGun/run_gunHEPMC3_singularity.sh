#!/bin/bash

source /cvmfs/cms.cern.ch/cmsset_default.sh
cmssw-el7 -- "source /cvmfs/sw.hsf.org/spackages6/key4hep-stack/2022-12-23/x86_64-centos7-gcc11.2.0-opt/ll3gi/setup.sh && ./gunHEPMC3 "$1" "
