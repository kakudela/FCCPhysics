# Analyzing events for FCC-ee

In this tutorial, we're going to analyze events from FCC-ee and try to measure some cross-sections. We'll learn:

- Part I: Basic analysis and event selection, plotting and cross-section measurement
- Part II: Fitting of histograms
- Part III: Jet clustering and flavor tagging
- PArt IV: Machine learning using Boosted Descsion Trees


# Prerequisites

#### Computing environment

Terminal, login, etc.

#### FCCAnalysis framework
The framework can be obtained following the instructions described [here](https://github.com/HEP-FCC/FCCAnalyses):

	## Check with Juraj -- can be done via cvmfs?
	cd go/to/my/directory
	source /cvmfs/sw.hsf.org/key4hep/setup.sh -r 2024-03-10
	git clone --branch pre-edm4hep1 git@github.com:HEP-FCC/FCCAnalyses.git
	cd FCCAnalyses
	source ./setup.sh
	fccanalysis build -j 8

The steps above have to be executed once. If you open a new terminal, you just need to repeat

	cd go/to/my/directory/FCCAnalyses
	source setup.sh



#### CMS Combination Tool (Combine)
We'll use the CMS Combination Tool for performing the statistical fits. This tool is built on top of the RooFit framework and is widely used within CMS for signal extraction, limit setting, and uncertainty estimation.

To simplify setup, we've pre-compiled a standalone version that can be run inside a Singularity image:

	## To be copied to central area
	/eos/experiment/fcc/users/j/jaeyserm/software/singularity/combine-standalone_v9.2.1.sif # CERN
	/work/submit/jaeyserm/software/docker/combine-standalone_v9.2.1.sif # MIT

Instructions on how to run it with Singularity are provided below.

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

	fccanalysis run FCCPhysics/tutorials/ZmumuHbb/analysis_ZmumuH.py

This script produces a ROOT file for each process in the directory `output/tutorial/zh_xsec/`, containing the corresponding histograms. You can inspect these files to see the cuts that have been applied.


The processes that are included in during the analysis are specified in the processList dictionary:

    processList = {
        'wzp6_ee_mumuH_ecm240': {'fraction': 1},
        'p8_ee_ZZ_ecm240': {'fraction': 1},
        'p8_ee_WW_mumu_ecm240': {'fraction': 1},
        'wzp6_ee_mumu_ecm240': {'fraction': 1},
    }

If you'd like to test the functionality of the script without processing the full dataset, you can reduce the number of events by adjusting the `fraction` values in `zh_xsec.py`. A fraction of `1` means all events are processed. You can reduce this value (e.g., to `0.1`) to speed up the creation of histograms for testing purposes.


Normalization can be automated during the analysis. To enable this, modify the following parameters in `z_mumu_xsec.py`:

    doScale = True
    intLumi = 10.8e6  # Integrated luminosity at the Higgs threshold (10.8 ab-1)

Run the script, and for each process a ROOT file is created in the specified output directory, containing all the histograms.


### Step 2: Plotting Histograms
To visualize the histograms and the cutflow plot, run the following command:

    fccanalysis plots FCCPhysics/tutorials/ZmumuHbb/plots_ZmumuH.py

Make sure to specify the output directory inside the `plots.py` script.

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

    fccanalysis combine FCCPhysics/tutorials/ZmumuHbb/combine_ZmumuH.py

The output files are saved in:

    output/ZmumuH/combine

This directory contains:
- A **text datacard**, also printed to screen  
- A **ROOT file** with the input histograms

The datacard is the key input to Combine — it tells the tool how to interpret the histograms in the context of a likelihood fit.


## Running the Fit and Output
To perform the likelihood fit using the CMS Combine tool inside a Singularity container, run:

    singularity exec --bind /work:/work /work/submit/jaeyserm/software/docker/combine-standalone_v9.2.1.sif bash -c 'cd output/ZmumuH/combine; text2workspace.py datacard.txt -o ws.root; combine -M MultiDimFit -v 10 --rMin 0.9 --rMax 1.1 --setParameters r=1 ws.root'

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












# Part II: Higgs to a pair of b-quarks
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
            'wz3p6_ee_mumuH_Hinv_ecm240'
        ],
        'bkg': [
            'wz3p6_ee_mumu_ecm240',
            'wz3p6_ee_tautau_ecm240',
            'p8_ee_WW_mumu_ecm240',
            'p8_ee_ZZ_ecm240'
        ]
    }

This setup treats only H → bb̄ as signal. All other Higgs decays (e.g. H → cc̄, WW, gluons, taus, etc.) are treated as an additional background (ZHnoBB), along with non-Higgs processes. 

Let's run the Combine fit using this new signal/background definition (change the process definitions in ```FCCPhysics/tutorials/ZmumuHbb/combine_ZmumuH.py```). Once you extract the uncertainty on the cross section for H → bb̄:

- What value do you obtain?
- Compare it to the previous inclusive Higgs result (0.6% uncertainty).
- Is the result compatible with the expected statistical worsening by a factor of 1/√0.58?


While we might expect a statistical degradation of the uncertainty by a factor of 1/√0.58 ~ 1.31 (this would naively lead to ~0.6% / 0.78 ≈ 0.8% uncertainty), in practice, the result may be worse. Why?

Because the non-bb Higgs decays form a significant background to H → bb̄ in our current selection. In order to overcome this, we will refine our selection and focus exclusively on the bb decays.





Both b-quarks from the Higgs decay manifest as jets in the detector — collimated sprays of particles resulting from hadronization. These final-state particles must be clustered into jets using a jet clustering algorithm. We use the [FastJet](https://indico.cern.ch/event/264054/contributions/592237/attachments/467910/648313/fastjet-doc-3.0.3.pdf) library, which provides a variety of jet algorithms commonly used in collider physics. For this analysis, we apply the Durham k<sub>T</sub> algorithm to cluster each event into exactly two jets, corresponding to the two b-quarks from the Higgs decay.

To distinguish b-jets from other types of jets (e.g. from gluons, c-quarks, or taus), we rely on jet flavor tagging. This process uses the properties of the jet — such as displaced vertices, track multiplicity, and invariant mass — to assign a probability that a given jet originated from a b-quark. This is done using a machine learning–based flavor tagger, specifically optimized for the FCC-ee environment. You can read more about the FCC-ee flavor tagger here: [FCC-ee Flavor Tagger – arXiv:2202.03285](https://arxiv.org/abs/2202.03285).


Both the jet clustering and flavor tagging can be easily performed within the FCCAnalyzer framework. To build on the previous `ZmumuH` analysis, we use a modified script:

    FCCPhysics/tutorials/ZmumuHbb/analysis_ZmumuHbb.py

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

    fccanalysis run FCCPhysics/tutorials/ZmumuHbb/combine_ZmumuHbb.py
    fccanalysis plots FCCPhysics/tutorials/ZmumuHbb/combine_ZmumuHbb.py
    fccanalysis combine FCCPhysics/tutorials/ZmumuHbb/combine_ZmumuHbb.py

Make sure to adapt the output directory in the script if needed. After running the full chain:

- How much does the significance improve after applying the b-tagging probability cut?
- How does the fit result (uncertainty on the H → bb̄ cross section) compare to the previous result without tagging?
- Is it closer to the expected statistical limit?
- Is the background from non-bb Higgs decays better suppressed? What about the backgrounds, in particular WW?

Use these comparisons to understand how flavor tagging enhances the precision of the measurement.








# Part III: Boosted Descision Tree (in progress)

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