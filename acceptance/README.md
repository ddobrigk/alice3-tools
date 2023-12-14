# Acceptance study tools

This directory contains a few files to aid in calculating the acceptance of strange decays in multiple tracker layouts. The files include: 

* `runacceptancetest.C`: a macro to be run compiled that allows for the calculation of baseline efficiencies given an arbitrary array of tracker layer radii. This macro is called by `runall.sh` using multiple magnetic field configurations. 
* `runall.sh`: a script to run acceptance checks at 0.5, 1.0, 1.5 and 2.0T by invoking `runacceptancetest.C` multiple times. Meant to be invoked after changing tracker layer radii in the `runacceptancetest.C` macro. 
* `DrawEfficiencies01.C`: a sample drawing macro to draw a few efficiencies and their ratios with respect to some baseline.
* `kinegen`: this directory contains the base code for generating kinematics trees for multiple magnetic field layouts.
