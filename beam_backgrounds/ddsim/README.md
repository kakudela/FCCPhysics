# Detector simulation with ddsim on the batch


We will use HTCondor to submit our jobs. More info on batch computing at MIT: https://submit.mit.edu/submit-users-guide/running.html#

Storage directory at MIT: 
    
    /ceph/submit/data/user/<first letter>/<USER>
    /ceph/submit/data/user/j/jaeyserm

To run over e.g. guinea pig files:

    python submit.py --gp_dir /ceph/submit/data/group/fcc/ee/detector/guineapig/FCCee_Z_4IP_04may23_FCCee_Z128/ --geometry_file $K4GEO/FCCee/CLD/compact/CLD_o2_v05/CLD_o2_v05.xml --storagedir /ceph/submit/data/user/j/jaeyserm/beam_backgrounds/


Note that the output will be stored in a directory with the name of the XML file, e.g.

    /ceph/submit/data/user/j/jaeyserm/beam_backgrounds/CLD_o2_v05/

So if you run on different detectors, please change the geometry file name!

