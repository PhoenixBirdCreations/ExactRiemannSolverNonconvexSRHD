# Exact Riemann solver for nonconvex SRHD with GGL EoS

This folder contains an implementation in C language of an exact Riemann solver for nonconvex Special Relativistic Hydrodynamics (SRHD) with a Gaussian-Gamma Law (GGL) equation of state (EoS).

Copyright 2023, Marina Berbel, All rights reserved

The code implements the solver in the paper 
## Exact Riemann solver for nonconvex relativistic hydrodynamics; Berbel, Serna, Marquina; JFM 2023

This folder contains:
### Libraries
Functions implemented for the solver, separated by topic for readability.

* auxFun: auxiliar functions for basic math
* EOS_GGL: all functions related to the EoS 
* hugoniot_curves: functions of Hugoniot curves
* integral_curves: functions of integral curves
* mixed_curves: functions of mixed curves
* principalLoops: functions monitoring the extremes of the curves, which should be calculated next and their intersection
* read_eos_par: reads the EoS par file
* read_ic_par: read the initial condition par file

### Main execution script
Reads the par files and start the solver. Outputs a log of the process in the console. Output the wave curves and the exact solution in txt files.

### Use of the code 
Download all files and put them in the same directory.

There is a Makefile distributed with the source code, simply change the name of the C compiler, if you use a different one, and type make.
It needs no external libraries besides the standard math, stdlib, string and stdio.

The executable is named exactRP and runs without any parameters. The details of the EoS are specified in a par file, eos.par. An example is provided. Two more examples, with GGL1 and GGL2 are included.

The initial conditions are specified in an ic.par file. An example is provided. Four more examples of initial conditions are included.

Once the problem is solved, the program creates two files with the names specified to solve the exact solution and the wave curves calculated.

The wave curves are stored in columns in the format

1-velocity &nbsp; &nbsp; &nbsp; 2-preassure &nbsp; &nbsp; &nbsp; 3-density &nbsp; &nbsp; &nbsp;   4-wave speed &nbsp; &nbsp; &nbsp; 5-characteristic speed &nbsp; &nbsp; &nbsp; 6-nonlinearityfactor

Different wave curves are separated by two line jumps. Heading each curve there is the text "L-X" or "R-X", for waves to the left ( L ) and the right ( R ), where X is H, I or M describing the wave type.

The exact solution is stored in columns in the format

1-position &nbsp; &nbsp; &nbsp; 2-velocity &nbsp; &nbsp; &nbsp; 3-preassure &nbsp; &nbsp; &nbsp; 4-density

First there are the waves to the left, from x=0 to the contact discontinuity. Then there are two line jumps and the waves to the right, from x=1 to the contact discontinuity. Finally after two line jumps there are the two states of the sides of the contact discontinuity.

### FAQs

**Fail in the execution due to v>1**

If the goal values for pressure and density are very far from the actual solution it is possible to obtain v>1 during the calculation of the waves due to floating point arithmetic errors. To solve this issue, it is enough to relax the goal values, as the code will modify them if necessary.

**Fail in mixed curves**

The termination state of a mixed curve that ends as a sonic shock is difficult to calculate since the equation to solve becomes very stiff. The rootfinder is calibrated for the examples proposed, but it is possible that for a given Riemann problem the solver fails. Please write to me at marinaberbelp@gmail.com and we can try to find a solution.



