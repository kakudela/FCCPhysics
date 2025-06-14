# Detector simulation with ddsim on the batch


We use **HTCondor** to submit our jobs. For more information on batch computing at MIT, refer to:  
[submit.mit.edu user guide](https://submit.mit.edu/submit-users-guide/running.html#)

Storage directory at MIT: 

    /ceph/submit/data/user/<first letter>/<USER>
    /ceph/submit/data/user/j/jaeyserm

To run over e.g. guinea pig files:

    python submit.py --gp_dir /ceph/submit/data/group/fcc/ee/detector/guineapig/FCCee_Z_4IP_04may23_FCCee_Z128/ --geometry_file $K4GEO/FCCee/CLD/compact/CLD_o2_v05/CLD_o2_v05.xml --storagedir /ceph/submit/data/user/j/jaeyserm/beam_backgrounds/


Note that the output will be stored in a directory with the name of the XML file, e.g.

    /ceph/submit/data/user/j/jaeyserm/beam_backgrounds/CLD_o2_v05/

So if you run on different detectors, please change the geometry file name!



# Changing geometry files

The geometry of the detector is defined by several XML files, together with C++ code that describes the full detector setup, including materials and components. For example, the CLD detector geometry is specified in the following file:

     $K4GEO/FCCee/CLD/compact/CLD_o2_v05/CLD_o2_v05.xml

Here, `$K4GEO` is an environment variable set when you source the `Key4HEP` stack, and `CLD_o2_v05` refers to the version of the CLD geometry. Inside this directory, many sub-XML files describe the individual components of the detector in more detail.


To run with a custom geometry, you typically modify the XML files. For instance, let’s change the magnetic field of the CLD detector from 2T to 3T. First, create a local directory to store the modified geometry:

    mkdir CLD_o2_v05_3T

Next, copy all the original XML files into this new directory:

    cp $K4GEO/FCCee/CLD/compact/CLD_o2_v05/* CLD_o2_v05_3T/

Rename the main XML file to reflect the new magnetic field:

    mv CLD_o2_v05_3T/CLD_o2_v05.xml CLD_o2_v05_3T/CLD_o2_v05_3T.xml

Then, edit the magnetic field value in the XML file. Locate and change the line:

    constant name="SolenoidField" value="3*tesla"/>

to

    <constant name="SolenoidField" value="3*tesla"/>

To run a physics simulation using this modified geometry, execute the following:

    python submit.py --gp_dir /ceph/submit/data/group/fcc/ee/detector/guineapig/FCCee_Z_4IP_04may23_FCCee_Z128/ --geometry_file CLD_o2_v05_3T/CLD_o2_v05_3T.xml --storagedir /ceph/submit/data/user/j/jaeyserm/beam_backgrounds/

The output will be stored in:

    /ceph/submit/data/group/fcc/ee/detector/ddsim//CLD_o2_v05_3T/FCCee_Z_4IP_04may23_FCCee_Z256/

This path automatically reflects the modified geometry name `CLD_o2_v05_3T`.