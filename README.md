# Contents

This repository contains code that goes with the paper *System size and boundaries determine the patterning dynamics of attracting active particles* (Rombouts et al.)

If you have questions you are welcome to send them to Jan Rombouts (jan.rombouts@ulb.be). 

## Python module files

There are three Python files that contain functions that are reused throughout. One for the simulations, one for bifurcation diagrams and one with general utilities. 

## Notebooks

There are demo notebooks to reproduce the main numerical and computational results of the paper. This includes
- A demo notebook to perform a single numerical simulation of the system.
- A notebook that shows how to compute a one-parameter bifurcation diagram. 
- There is a separate notebook to visualize and browse bifurcation diagrams. Here, a saved bifurcation diagram (pickled file) can be loaded and then visualized.
- A notebook to browse through the experimental data (the density profiles). 
- A notebook that shows how we performed the pattern classification used in Fig. 5. 
- A notebook that shows how we fitted the experimental data in the initial regime to obtain estimates for $\bar\alpha$, which is used in Fig. 5 in the paper. 
- A notebook that demonstrates the main calculations in relation to the linear stability and the calculation of the eigenvalues (Fig. 2)

## Figure scripts

- All python scripts that generate the figures of the paper are included. 
- The data necessary to produce these figures is included in a separate folder. 
- there is a *utils.py* file in the figure scripts folder, which defines a few functions that are used in different figures, such as loading experimental profiles etc. 

## Data

There are two data folders, one is just called *data*, and contains all experimental density profiles. These are the processed profiles, the raw image files are shared with our other manuscript (in preparation). The other data folder is called *data_for_figures*, and contains all files that are necessary to produce the figures.
