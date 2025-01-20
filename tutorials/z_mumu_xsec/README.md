# Forward-Backward Asymmetry at the Z Pole
## Introduction
The vector and axial-vector components of the Z boson result in an asymmetric angular distribution of the outgoing leptons with respect to the polar angle, $\theta$. This asymmetry is characterized by a parameter known as $A_{FB}$ (forward-backward asymmetry), which can be quantified and measured. Achieving high-precision measurements of $A_{FB}$ is crucial, as it provides valuable insights into the underlying physics. Therefore, the value of $A_{FB}$ can be slightly different as each event generator has its own implementation of the physics model, as we will see in this tutorial. Notably, $A_{FB}(s)$ typically varies with the center-of-mass energy $s$.

The parameterization of the forward-backward asymmetry solely depends on the angles of both leptons with respect to the incoming particles (theta). One defines the angle $\cos(\theta_c)$ as follows:

$$
\cos\theta_c = \frac{\sin{(\theta_+-\theta_-)}}{\sin{\theta_+}+\sin{\theta_-}},
$$
where $\theta_+$ and $\theta_-$ represent the polar angles of the positive and negative muons, respectively, measured in the lab frame. The angle $\theta_c$ is the scattering angle of the negative muon in the reduced center-of-mass frame of the muon pair, assuming that the initial state radiated photons have zero transverse momenta. The differential cross-section for muon pair production is given by:

$$
\frac{d\sigma(s)}{d\cos\theta_c} = \sigma(s) \times \left( \frac{3}{8}(1+\cos^2\theta_c) + A_{FB}(s) \cos\theta_c \right).
$$
The asymmetry is defined as the asymmetry of the $\cos\theta_c$ distribution:
$$
A_{FB} = \frac{\sigma_F - \sigma_B}{\sigma_F + \sigma_B} = \frac{N_F - N_B}{N_F + N_B},
$$
with $N_F$ and $N_B$ the number of events counted in the forward and backward regions, respectively (corresponding to the number of events in the positive and negative sides of the $\cos\theta_c$ distribution).

## Analysis
In this tutorial, we extract the forward-backward asymmetry at the Z pole in di-muon events for different event generators (Whizard, KKKMC, and Pythia). The analysis consists of the following tasks:

1. Implement the angle $\cos\theta_c$ in the analyzer and plot it for the different generators. Is the plot symmetric?
2. Calculate $A_{FB}$ by integrating over the positive and negative $\cos\theta_c$ distribution for the different generators.
3. Derive the statistical uncertainty of $A_{FB}$.
4. Rather than counting the numbers on both sides of the spectrum, perform a $\chi^2$ fit on the $\cos\theta_c$ distribution according to the formula of the differential cross-section above, in order to extract $A_{FB}$. Does it agree with the $A_{FB}$ using the integration method above? Interpret with the associated uncertainties and compare with LEP and FCC-ee integrated luminosities.
5. Compare the different event generators and their value of $A_{FB}$.
6. Extract the gen-level $A_{FB}$ and compare with the reco-level $A_{FB}$.

To run the forward-backward asymmetry analysis, execute the following script from the main `FCCAnalyzer` directory (to quickly run over a few files, adjust the `fraction` in the `z_mumu_afb.py` file):

```shell
fccanalysis run FCCPhysics/tutorials/z_mumu_afb/z_mumu_afb.py
```

This produces a ROOT file per process in the directory `output/tutorials/z_mumu_afb/`, which contains the histograms. To plot and fit the forward-backward asymmetry, run the following command in Jupyter Notebook:

```shell
FCCPhysics/tutorials/afb/afb_fit.ipynb
```

Alternatively, a standalone script written in ROOT is available to extract the forward-backward asymmetry:

```shell
python FCCPhysics/tutorials/afb/afb_fit.py --proc wzp6_ee_mumu_ecm91p2
python FCCPhysics/tutorials/afb/afb_fit.py --proc kkmcee_ee_mumu_ecm91p2
python FCCPhysics/tutorials/afb/afb_fit.py --proc p8_ee_Zmumu_ecm91
```
