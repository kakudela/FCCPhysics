# Advanced Techniques: ROOT DataFrames

In this tutorial, you will learn how to use DataFrames to speed up your analysis code.
We cover ROOT's version of DataFrames, called RDFs.

# DataFrames, in general

As we talked about in the previous tutorial, our data is organized in separate events. We loop over these events, as we need to analyze each one individually, and then aggregate the results.

From a computing point of view, this is a achieved using `for` loops over the events. While being easy to write and understand, writing `for` loops in `PyROOT` is not a very efficient way to process the large datasets that we have to work with in this field.

Instead, we can use DataFrames, which are highly optimized frameworks to manipulate our data. The key point here is that you define actions that are applied on the whole data at once; for example, if you apply a selection, it will be performed on all events in the DataFrame. This type of columnar operation is a lot quicker computationally. (Don't be afraid, these are just `for` loops under the hood! They are just more efficient than the `for` loops in `PyROOT`, and take care of some parallelization for you, allowing for things to scale seamlessly, but more on that later.)

Here is a simple comparison of `for` loops and DataFrames approaches:

```python
# simple for loop
for event in events:
    for particle in event.particles:
        if particle.pt > 10:
            hist.fill(particle.pt)
```

```python
# RDF
df = df.Filter("particles.pt > 10")
df.Histo1D("particles.pt")
```

# ROOT DataFrames

ROOT's [docs](https://root.cern/doc/master/classROOT_1_1RDataFrame.html) provide good documentation, tutorials, and a crash course on RDFs, check them out!

For complicated operations, you will need to write functions that act on the data. This needs to be done in C++, because of how ROOT compiles the code. You can find some examples in the `functions.h` script of this repo. Further tutorials will provide more real-life examples of such functions. The `FCCAnalyses` framework already provides a lot of useful functions, so check out their documentation.

### Example

Here is a simple script that reads a reads a ROOT file, selects on MCParticles with at least 10 GeV of energy, and fills a histogram with the energies of these particles. This is the same analysis we performed in the previous tutorial, with a `for` loop, now done with an RDF.

In order to run this, create a file called `energy.h` with the following content:

```cpp
#include "FCCAnalyses/MCParticle.h"

ROOT::VecOps::RVec<float> get_energy(const ROOT::VecOps::RVec<edm4hep::MCParticleData>& particles) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: particles) {
    TLorentzVector tlv;
    tlv.SetXYZM(p.momentum.x, p.momentum.y, p.momentum.z, p.mass);
    float energy = tlv.E();
    if (energy > 10) {
      result.push_back(energy);
    }
  }
  return result;
}
```

This function calculates the energy of the particles, and returns a vector with the energies of the particles with more than 10 GeV.
This function is then used in the RDF to fill the histogram.

```python
import ROOT

ROOT.gInterpreter.ProcessLine("""
#include "energy.h"
"""
)

# Open the file
f = ROOT.TFile.Open("/ceph/submit/data/group/cms/store/fccee/samples/winter2023/CLD_FullSim/wzp6_ee_mumuH_ecm240/mumuH_rec_16610_191.root")

# Get the tree
tree = f.Get("events")

# Create the RDF
df = ROOT.RDataFrame(tree)

# Define the energy of the particles, and save if above 10 GeV
df = df.Define("Energy", "get_energy(MCParticles)")

# Fill the histogram with the calculated energies
# We defined 100 bins between 0 and 1000 GeV
h = df.Histo1D(("MCParticles_Energies", "", *(100, 0, 1000)), "Energy")

# Draw the histogram
c = ROOT.TCanvas()
h.Draw()
c.Draw()

# save the histogram
c.SaveAs("energy.png")
```

> *Exercise*: run this script, and produce the energy histogram. How does it compare with the one from the previous tutorial? Write the analysis to extract the invariant mass of the dimuon system using RDFs. Try running over more events, and note the speed between the simple `for` loop and RDFs.