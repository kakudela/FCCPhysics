___________________________________________________________________________________
___________________________________________________________________________________

GUINEA-PIG TIMING STUDIES

___________________________________________________________________________________
___________________________________________________________________________________


env.sh: source

acc.dat: contains all of the beam parameters and Guinea-Pig configs


___________________________________________________________________________________
___________________________________________________________________________________


run_gp_timing.sh

Automatically runs different Tracking/noTracking and gridsN of Guinea-Pig for
ease because I'm lazy, writes timing information in a csv that is analyzed and
plotted by another script.


___________________________________________________________________________________
___________________________________________________________________________________


plot_gp_timings.py

Takes the CSV that was written after the timing studies, plots number of GP grids
vs. computing time

- Also makes some plots based on the grid area/volume configurations, rates etc.
- It turns out that all of the cells in each grid actually do run the computations,
  even if, eg., the volume in grid five that overlaps with grid 4 is never seen by IPCs
    - This is due to the Fast Fourier Transform package; all cells need to be
      computed


___________________________________________________________________________________
___________________________________________________________________________________
