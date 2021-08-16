# NRL MAG NOH PROBLEM 

This is supplementary code for the JFM paper "Stable and unstable supersonic stagnation of an axisymmetric rotating magnetized plasma"

Requirements:
1) for generating self-similar solutions - Python 3 (preferably with numba)
2) for verification with MHD code - Athena++

To reproduce figures from the paper, run:

Fig 2 bottom

python3 ./verify_new.py 1.5 1.66667 1e5 1e5 6e-5

Fig 2 top

python3 ./verify_new.py 2.6 1.1 1.198e3 1.198e3 3.256e-9

Fig 6:

python3 ./verify_new.py 2.0 1.66667 1e5 1e5 6e-5 --rot=0.02

Fig 5, bottom:

python3 ./verify_new.py 1.5 1.66667 1e5 3e4 1e-4 --v0solve --nr=10000

Last fig, asympt:

python3 ./asympt.py 1.5 1.2 1 1 0.6 --u0max=100

Initial pressure example:

python3 ./verify_new.py 1.5 1.66667 1e5 1e5 6e-5 --beta=20
