# Higgs Total Cross-Section

In this tutorial, we extract the total Higgs cross-section using the recoil method. A cutflow analysis is set up to reduce the main backgrounds (WW, ZZ, and Z/γ). The final recoil histogram is then used with the CMS Combine tool to extract the uncertainty on the cross-section through a binned maximum likelihood fit.

## Step 1: Event Selection and Histogram Creation

To generate the event selection and histograms for all signal and background samples, run the following script from the main FCCAnalyses directory:

```shell
fccanalysis run FCCPhysics/tutorials/zh_xsec/zh_xsec.py
```

This script produces a ROOT file per process in the directory `output/tutorial/zh_xsec/`, containing the histograms.

If you want to quickly test the functionality of the script without running over the full set of events, you can adjust the fraction of files in the `zh_xsec.py` file:

```shell
processList = {
    'wzp6_ee_mumuH_ecm240': {'fraction':1},
    'p8_ee_ZZ_ecm240': {'fraction':1},
    'p8_ee_WW_mumu_ecm240': {'fraction':1},
    'wzp6_ee_mumu_ecm240': {'fraction':1},
}
```

A fraction of `1` means all events are processed.

## Step 2: Plotting Histograms

To visualize the main histograms and the cutflow histogram, execute the following command:

```shell
fccanalysis plots FCCPhysics/tutorials/zh_xsec/plots.py
```

Make sure to specify the output directory in the `plots.py` script.

## Step 3: Preparing for the Fit

To extract the uncertainty on the cross-section, we first need to prepare the datacards. Use the following command to generate the necessary text and ROOT files with the final histograms for the fit:

```shell
fccanalysis combine FCCPhysics/tutorials/zh_xsec/combine.py
```

The output files are saved in the directory `output/tutorial/zh_xsec/combine`.

## Step 4: Running the Fit with CMS Combine

The next step is to perform the fit using the CMS Combine tool. Run the following command in a Singularity container:

```shell
singularity exec --bind /work:/work /work/submit/jaeyserm/docker/combine-standalone_v9.2.1.sif bash -c 'cd output/tutorial/zh_xsec/combine; text2workspace.py datacard.txt -o ws.root; combine -M MultiDimFit -v 10 --rMin 0.9 --rMax 1.1 --setParameters r=1 ws.root'
```

## Step 5: Extracting the Results

After running the fit, the output will look similar to the following:

```shell
Minuit2Minimizer : Valid minimum - status = 0
FVAL  = 0
Edm   = 2.38118773825041254e-12
Nfcn  = 27
bkg_norm          = 0    +/-  0.0957588
r         = 1    +/-  0.00755054        (limited)
Minimization finished with status=0
Minimization success! status=0
Minimized in 0.036844 seconds (0.030000 CPU time)
FINAL NLL - NLL0 VALUE = 6.982164185e-07
```

From this output, you can extract the parameter r and its uncertainty. In this example, the result is r = 1 ± 0.00755054, corresponding to a relative uncertainty of 0.755%.