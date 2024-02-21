# NRL MAG NOH PROBLEM 

NRL Mag Noh is a five-parametric family of self-similar solutions of ideal MHD. These exact solutions which describe implosion in cylindrical geometry were originally intended to test MHD codes, but can also be used to study basic physics of implosion flows, stagnation, shocks, shock stability, magnetorotational instability and some. Some aspects of NRL Mag Noh were described in the Journal of Fluid Dynamics paper "Stable and unstable supersonic stagnation of an axisymmetric rotating magnetized plasma" by Andrey Beresnyak, Alexander L. Velikovich, John L. Giuliani, and Arati Dasgupta: https://doi.org/10.1017/jfm.2022.77.

The infamous *Problem 3* that broke almost every ideal MHD code in Cartesian geometry can be reproduced by running:

`verify.py 2.6 1.1 1.198e3 1.198e3 3.256e-9`

This generates a tabulated self-similar solution, as well as runs Athena++ and produces the verificaton figure.

Below are the command line options and software requirements to generate the JFM paper figures using our python scripts.

Requirements:
1) for generating self-similar solutions - Python (preferably with numba)
2) for verification with MHD code - Athena++

To compile Athena for verification purposes, use magnoh.cpp from this directory. The currently publicly released magnoh.cpp from Athena pgen directory does not set up rotation.

To compile Athena for stability analysis, use magnoh_cyl2.cpp

To see the options use --help or -h option of python scripts

To reproduce figures from the paper, have "athena" executable in the same directory as the python script and run:

`mkdir solutions`

Fig 2 bottom:
`verify.py 1.5 1.66667 1e5 1e5 6e-5`

Fig 2 top:
`verify.py 2.6 1.1 1.198e3 1.198e3 3.256e-9`

Fig 6:
`verify.py 2.0 1.66667 1e5 1e5 6e-5 --rot=0.02`

Fig 5, bottom:
`verify.py 1.5 1.66667 1e5 3e4 1e-4 --v0solve --nr=10000`

Last fig, asympt:
`asympt.py 1.5 1.2 1 1 0.6 --u0max=100`

Initial pressure example:
`verify.py 1.5 1.66667 1e5 1e5 6e-5 --beta=20`

The figures and the text file with tabulated solutions will be generated in the "solutions" subdirectory.

Rotating perturbed solution:
![rotating perturbed](https://github.com/beresnyak/magnoh/blob/main/examples/rotating_perturbed1.png)



Layman's introduction of what 'self-similar' means:
![self-similar plot](https://github.com/beresnyak/magnoh/blob/main/examples/ss_plot.png)
