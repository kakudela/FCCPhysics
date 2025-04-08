
# Analyzing events for FCC-ee
In this tutorial, we're going to analyze events from the **FCC-ee** and measure cross-sections for selected processes. Throughout the tutorial, you'll learn how to perform a basic analysis, fit histograms, apply jet clustering and flavor tagging, and use machine learning techniques:

- **Part I:** Basic analysis and event selection, plotting, and cross-section measurement  
- **Part II:** Fitting of histograms  
- **Part III:** Jet clustering and flavor tagging  
- **Part IV:** Machine learning using Boosted Decision Trees



# Prerequisites

#### Computing environment

Instructions for account creation, logging in into the SubMIT cluster can be found [here](https://github.com/jeyserma/FCCPhysics/blob/main/tutorials/ZmumuHbb/SUBMIT.md).



#### FCCAnalysis framework
We will be working with the **FCCAnalyses** framework, available on [GitHub](https://github.com/HEP-FCC/FCCAnalyses). This is a common analysis framework developed for FCC-related studies. It lets you run full analyses over existing simulated samples, apply event selections, and produce plots and histograms. The input samples used in this tutorial are available [here](https://fcc-physics-events.web.cern.ch/fcc-ee/delphes/winter2023/idea/).


If you are running the tutorial on the MIT computing infrastructure, there is a pre-installed version of the analysis software available. You can set it up by running (this command must be executed every time you log into a new terminal session):

    source /work/submit/jaeyserm/software/FCCAnalyses/setup.sh





If you're running at CERN or prefer to install the analysis framework yourself, follow the instructions below (adapted from the [FCCAnalyses GitHub repository](https://github.com/HEP-FCC/FCCAnalyses)):

    cd go/to/my/directory
    source /cvmfs/sw.hsf.org/key4hep/setup.sh -r 2024-03-10
    git clone --branch pre-edm4hep1 git@github.com:HEP-FCC/FCCAnalyses.git
    cd FCCAnalyses
    source ./setup.sh
    fccanalysis build -j 8

These setup steps need to be executed only **once** during the initial installation.

For subsequent sessions, simply run:

    cd go/to/my/directory/FCCAnalyses
    source setup.sh

####  Getting the Tutorial Files

We’ve prepared a few Python files to guide you through this tutorial. You can find them in a separate repository  
[here](https://github.com/jeyserma/FCCPhysics/tree/main/tutorials/ZmumuHbb) [to be defined where we will store the tutorial files]

You can either download the files individually from the repository or clone the entire repository using:

    git clone https://github.com/jeyserma/FCCPhysics.git
    cd tutorials/ZmumuHbb




#### CMS Combination Tool (Combine)
We'll use the CMS Combination Tool for performing the statistical fits. This tool is built on top of the RooFit framework and is widely used within CMS for signal extraction, limit setting, and uncertainty estimation.

To simplify setup, we'll use a pre-compiled standalone version that can be run inside a Singularity image — an isolated Linux environment where we can install and run software without affecting the host system. This helps ensure the analysis runs the same way on different machines.

You can find the pre-built image at the following locations:

	/work/submit/jaeyserm/software/docker/combine-standalone_v9.2.1.sif # MIT
	/eos/project/f/fccsw-web/www/analysis/auxiliary/combine-standalone_v9.2.1.sif # CERN  

Instructions on how to run it with Singularity are provided below in the tutorial.

More information on how to compile the package locally and use all its features can be found in the [official documentation](https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/latest/#oustide-of-cmssw-recommended-for-non-cms-users).




# Part I: Measuring a Cross-Section
In this first part of the tutorial, we will perform a basic event selection at the Higgs threshold. Our target process is:

$$
	\rm e⁺+e⁻ → ZH
$$
with the Z boson decaying into a pair of muons (Z → μ⁺μ⁻), and the Higgs boson allowed to decay into any final state.

We aim to extract the total Higgs production cross-section using the recoil method. To isolate the signal, we perform a cutflow analysis to reduce the main backgrounds, which include:

- WW production  
- ZZ production  
- Z/γ* (Drell-Yan-like processes)

Finally, we construct a recoil mass histogram, which is used as input to the [CMS Combine tool](https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/) to estimate the uncertainty on the cross-section via a binned maximum likelihood fit.


### Step 1: Running the Event Selection and Producing Histograms

To generate the event selection and histograms for all signal and background samples (which are pre-generated using the FastSim IDEA detector), run the following script from the main `FCCAnalyses` directory:

	fccanalysis run analysis_ZmumuH.py

You can inspect this file to see the cuts that have been applied. The processes that are included in during the analysis are specified in the processList dictionary:

    processList = {
        'wzp6_ee_mumuH_HXX_ecm240':      {'fraction': 1},
        'p8_ee_ZZ_ecm240':               {'fraction': 1},
        'p8_ee_WW_mumu_ecm240':          {'fraction': 1},
        'wzp6_ee_mumu_ecm240':           {'fraction': 1},
    }

For the Higgs signal sample `wzp6_ee_mumuH_HXX_ecm240`, we have several samples available for each individual Higgs decay .

If you'd like to test the functionality of the script without processing the full datasets, you can reduce the number of events by adjusting the `fraction` values in `zh_xsec.py`. A fraction of `1` means all events are processed. You can reduce this value (e.g., to `0.1`) to speed up the creation of histograms for testing purposes.


Normalization can be automated during the analysis. To enable this, modify the following parameters in `z_mumu_xsec.py`:

    doScale = True
    intLumi = 10.8e6  # Integrated luminosity at the Higgs threshold (10.8 ab-1)

After running the script, a ROOT file is created for each process in the output directory (`output/ZmumuH/histmaker/`). These files contain all the histograms needed for the next steps of the analysis. If you run into issues executing the setup or want to skip the processing step, we’ve pre-generated the files with full statistics. You can access them here:

	/ceph/submit/data/group/fcc/ee/tutorials/FNAL2025/ZmumuH/histmaker/ # MIT
	/eos/project/f/fccsw-web/www/tutorials/FNAL2025/ZmumuH/histmaker/ # CERN




### Step 2: Plotting Histograms
To visualize the histograms and the cutflow plot, run the following command:

    fccanalysis plots plots_ZmumuH.py

Make sure to specify the input and output directories inside the `plots.py` script (e.g. if you want to plot and run the pre-generated files as explained above, you should point to that directory).

The **cutflow plot** provides a clear illustration of how background events are progressively reduced by each selection cut, while ideally retaining the signal. To quantify the efficiency of this selection, we calculate the significance:
$$
\rm significance = \frac{S}{\sqrt{S + B}},
$$

where `S` and `B` are the number of signal and background events, respectively, after each cut. As background is reduced and signal is preserved, the significance increases. This value reflects how confidently we can measure the signal in the presence of background. The uncertainty, in %, on the signal is calculated as 1/significance. In the special case where background is negligible, this expression simplifies to:
$$
\rm significance = \sqrt{S},
$$

and the **uncertainty** on the signal yield becomes:
$$
\rm uncertainty =\frac{1}{\sqrt{S}}.
$$

Maximizing the significance and mimnimizing the uncertainty is a key goal in any physics analysis — it's our job as physicists! More advanced methods, such as Machine Learning (covered later), can help you push this even further.

The significance after each cut has already been calculated for you and is listed in the file `cutFlow.txt`. Do you observe a final significance of 120? What is the corresponding uncertainty on the signal process?

Now try tightening the recoil mass cut to see if you can further improve the significance. How much more can you squeeze out of it?




# Part II: Statistical Analysis
From the definition of significance above, we can also evaluate it **bin-by-bin**, treating each histogram bin as a separate measurement. Instead of tightening the recoil mass cut, we can instead **combine information from all bins** to extract a **single uncertainty** on the signal yield.

This is achieved using a **likelihood fit**, which does this combination optimally. The advantage of likelihood fits is that they also allow you to assign **systematic uncertainties** on both signal and background components — although that’s beyond the scope of this tutorial.

## Preparing the Datacards
To extract the uncertainty on the signal cross-section, we first need to prepare the Combine-compatible datacards. Use the following command to generate the required text and ROOT files containing the final histograms:

    fccanalysis combine combine_ZmumuH.py

The output files are saved in:

    output/ZmumuH/combine

This directory contains:
- A **text datacard**, also printed to screen  
- A **ROOT file** with the input histograms

The datacard is the key input to Combine — it tells the tool how to interpret the histograms in the context of a likelihood fit.


## Running the Fit and Output
To perform the likelihood fit using the CMS Combine tool inside a Singularity container, run (MIT and CERN respectively):

    singularity exec --bind /work:/work /work/submit/jaeyserm/software/docker/combine-standalone_v9.2.1.sif bash -c 'cd output/ZmumuH/combine; text2workspace.py datacard.txt -o ws.root; combine -M MultiDimFit -v 10 --rMin 0.9 --rMax 1.1 --setParameters r=1 ws.root'
    singularity exec /eos/project/f/fccsw-web/www/analysis/auxiliary/combine-standalone_v9.2.1.sif bash -c 'cd output/ZmumuH/combine; text2workspace.py datacard.txt -o ws.root; combine -M MultiDimFit -v 10 --rMin 0.9 --rMax 1.1 --setParameters r=1 ws.root'



After the fit runs, you’ll see output similar to this:

    Minuit2Minimizer : Valid minimum - status = 0
    FVAL  = -2.77185177588457066e-10
    Edm   = 3.21698353423223593e-13
    Nfcn  = 15
    r         = 1    +/-  0.00607216        (limited)
    Minimization finished with status=0
    Minimization success! status=0
    Minimized in 0.043089 seconds (0.040000 CPU time)
    FINAL NLL - NLL0 VALUE = -2.771851776e-10

From this output, you can extract:
- The best-fit value of `r` (signal strength), which is `1.0`
- The uncertainty on `r`, which is approximately `±0.006`, or 0.6%


How Does This Compare to the Significance?
- The **significance** gives a measure of how likely it is that the observed signal stands out from the background — useful for discovery.
- The **likelihood fit** provides the uncertainty on the signal strength, which is key for precision measurements.

Is 0.6% uncertainty better or worse than what you got from the cut-based significance? Why?
Think about:
- Whether combining bins gives you more statistical power  
- How much information you lost by applying a single tight cut  
- The advantages of fitting a shape vs. counting events












# Part III: Higgs to a pair of b-quarks

## Extract the uncertainty for H → bb̄ 
The Higgs boson predominantly decays to a pair of b-quarks, with a branching ratio of approximately 58%. One of the key objectives of FCC-ee is to precisely measure the cross section for H → bb̄, which directly relates to the Higgs coupling to b-quarks.

So far, we have considered all Higgs decay modes (H → bb̄, cc̄, τ⁺τ⁻, WW, ZZ, etc.) as signal. To focus specifically on H → bb̄, we must redefine the signal** and background processes in the Combine setup:

    sig_procs = {'sig': ['wzp6_ee_mumuH_Hbb_ecm240']}
    bkg_procs = {
        'ZHnoBB': [
            'wzp6_ee_mumuH_Hcc_ecm240',
            'wzp6_ee_mumuH_Hss_ecm240',
            'wzp6_ee_mumuH_Hgg_ecm240',
            'wzp6_ee_mumuH_Haa_ecm240',
            'wzp6_ee_mumuH_HZa_ecm240',
            'wzp6_ee_mumuH_HWW_ecm240',
            'wzp6_ee_mumuH_HZZ_ecm240',
            'wzp6_ee_mumuH_Hmumu_ecm240',
            'wzp6_ee_mumuH_Htautau_ecm240',
        ],
        'bkg': [
            'wzp6_ee_mumu_ecm240',
            'wzp6_ee_tautau_ecm240',
            'p8_ee_WW_mumu_ecm240',
            'p8_ee_ZZ_ecm240'
        ]
    }

This setup treats only H → bb̄ as signal. All other Higgs decays (e.g. H → cc̄, WW, gluons, taus, etc.) are treated as an additional background (ZHnoBB), along with non-Higgs processes. 

Let's run the Combine fit using this new signal/background definition: change the process definitions in ```combine_ZmumuH.py``` as defined above and re-run Combine. Extract the uncertainty on the cross section for H → bb̄:

- What value do you obtain?
- Compare it to the previous inclusive Higgs result (0.6% uncertainty).
- Is the result compatible with the expected statistical worsening by a factor of 1/√0.58?


While we might expect a statistical degradation of the uncertainty by a factor of 1/√0.58 ~ 1.31, that would naively lead to an 0.6%*1.31 ≈ 0.8% uncertainty), in practice, the result is worse. Why?

Because the non-bb Higgs decays form a significant background to H → bb̄ in our current selection, and they all exhibit the same shape of the recoil distribution. In order to overcome this, we will refine our selection and focus exclusively on the bb decays.



## Using the Flavor tagger
Both b-quarks from the Higgs decay manifest as jets in the detector — collimated sprays of particles resulting from hadronization. These final-state particles must be clustered into jets using a jet clustering algorithm. We use the [FastJet](https://indico.cern.ch/event/264054/contributions/592237/attachments/467910/648313/fastjet-doc-3.0.3.pdf) library, which provides a variety of jet algorithms commonly used in collider physics. For this analysis, we apply the Durham k<sub>T</sub> algorithm to cluster each event into exactly two jets, corresponding to the two b-quarks from the Higgs decay.

To distinguish b-jets from other types of jets (e.g. from gluons, c-quarks, or taus), we rely on jet flavor tagging. This process uses the properties of the jet — such as displaced vertices, track multiplicity, and invariant mass — to assign a probability that a given jet originated from a b-quark. This is done using a machine learning–based flavor tagger, specifically optimized for the FCC-ee environment. You can read more about the FCC-ee flavor tagger here: [FCC-ee Flavor Tagger – arXiv:2202.03285](https://arxiv.org/abs/2202.03285).


Both the jet clustering and flavor tagging can be easily performed within the FCCAnalyzer framework. To build on the previous `ZmumuH` analysis, we use a modified script:

    FCCPhysics run analysis_ZmumuHbb.py

Here are some important lines to inspect in comparison with the previous file:

    # Remove the muons from the event before clustering
    df = df.Define("rps_no_muons", "FCCAnalyses::ReconstructedParticle::remove(ReconstructedParticles, muons)")

    # Cluster the remaining particles into exactly 2 jets
    jetClusteringHelper = ExclusiveJetClusteringHelper("rps_no_muons", 2, "")
    df = jetClusteringHelper.define(df)  # Performs clustering and creates jet collections

    # Setup and run the flavor tagger
    df = jetFlavourHelper.define(df)
    df = jetFlavourHelper.inference(weaver_preproc, weaver_model, df)  # Run inference

The jet clustering produces a collection of reconstructed jets, which can be converted to Lorentz vectors for physics analysis (e.g. dijet invariant mass). The flavor tagger returns per-jet probabilities for the jet being a `b`, `c`, `s`, `g`, or `τ` jet.

To select b-jets, we apply the following cut:

    df = df.Filter("recojet_isB[0] > 0.5 && recojet_isB[1] > 0.5")

This keeps only events where both jets have a b-tagging probability greater than 0.5.

Execute the full analysis pipeline:

    fccanalysis run analysis_ZmumuHbb.py
    fccanalysis plots plots_ZmumuHbb.py
    fccanalysis combine combine_ZmumuHbb.py
    singularity exec --bind /work:/work /work/submit/jaeyserm/software/docker/combine-standalone_v9.2.1.sif bash -c 'cd output/ZmumuHbb/combine; text2workspace.py datacard.txt -o ws.root; combine -M MultiDimFit -v 10 --rMin 0.9 --rMax 1.1 --setParameters r=1 ws.root'
    singularity exec /eos/project/f/fccsw-web/www/analysis/auxiliary/combine-standalone_v9.2.1.sif bash -c 'cd output/ZmumuHbb/combine; text2workspace.py datacard.txt -o ws.root; combine -M MultiDimFit -v 10 --rMin 0.9 --rMax 1.1 --setParameters r=1 ws.root'

Make sure to adapt the input/output directory in the scripts if needed. 

As you will notice, the first command might take a while to run — this is because the tagger inference step is relatively slow. If you don't want to wait for its output, or if you run into issues executing the setup, you can skip this step entirely. We’ve pre-generated the files with full statistics, and you can access them here:

	/ceph/submit/data/group/fcc/ee/tutorials/FNAL2025/ZmumuHbb/histmaker/ # MIT
	/eos/project/f/fccsw-web/www/tutorials/FNAL2025/ZmumuHbb/histmaker/ # CERN

If you use these samples, make sure to update the input directory in the plotting and combine scripts.


After running the full chain:

- How much does the significance improve after applying the b-tagging probability cut?
- How does the fit result (uncertainty on the H → bb̄ cross section) compare to the previous result without tagging?
- Is it closer to the expected statistical limit?
- Is the background from non-bb Higgs decays better suppressed? What about the backgrounds, in particular WW?

Use these comparisons to understand how flavor tagging enhances the precision of the measurement.








# Part IV: Boosted Descision Tree (in progress)

Instread of cut and counting, we will use Machine Learning techniques to improve our result. We need to prepare variables tthat distinguish between signal and background. Candidates are:

- Recoil mass
- Higgs mass
- Jet kinematics
- Flavor information of the jets


### Preparing the input variables

### Run the Training and Evaluation

### Apply the inference

### Evaluate the performance

Evaluate the performance by fitting on the MVA discriminator. How does it does the result change w.r.t. the previous result?