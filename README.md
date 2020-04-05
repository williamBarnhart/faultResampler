# faultResampler
MATLAB software for inverting geodetic observations (InSAR, optical, GPS) for fault slip onto variably sized, triangular fault patches

FaultResampler is designed for fault-slip inversions from geodetic data (GPS, InSAR, optical offsets, etc.). The idea is simple: Inversion resolution (what you can actually say about fault slip) is fundamentally dependent on the quality of your data and the location of your data relative to the fault plane.  This idea plays the greatest part in where on the fault you can resolve slip well.  Near data points, you have good resolution of slip; far from data points, you have poor resolution.  FaultResampler iteratively discretizes an input fault model to produce a fault model and slip distribution that accurately reflect the resolution power of your data.

Reference: 
Barnhart, W.D., Lohman, R.B. (2010) Automated fault model discretization for inversions for coseismic slip distributions, J. Geophys. Res. 115, B10419.

## Requirements:

MATLAB,
MATLAB Optimization Toolbox (for running the functions lsqlin and lsqnonlin)

Brendan Meade's (Harvard University) triangular dislocation code: https://github.com/brendanjmeade/tde

## Installation:
1) Download the respository into a location where you store MATLAB scripts. Download the triangular dislocations code and place it in this folder as well.
2) Read the README.pdf file
3) Open MATLAB and add the program and subfolders to your path: addpath(genpath('/path/to/matlabfiles/faultResampler'))

Questions or bugs? Email Bill Barnhart (william-barnhart-1 at uiowa dot edu)
