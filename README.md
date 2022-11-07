# Lead-Green-s-functions-from-quadratic-eigenvalue-problems-without-mode-velocity-calculations
Accompanying code for the paper Lead Green's functions from quadratic eigenvalue problems without mode velocity calculations.

The folder \scr contains julia code that implements the methods described in the paper.The timing and forward error data used in the paper are in the
folder \Data.

It is possible to create similar data by running the Calculate_lead_energybands.ipynb and Create_data_for_Fig2.ipynb jupyter notebooks and the following julia scripts

```
julia --threads 1 --procs 1 Create_data_for_Fig3.jl
julia --threads 1 --procs 1 Create_data_for_Figs_5_and_6.jl
```

The folder \Scripts_to_create_figures contains python scripts to create the figures in the paper, stored in the folder \Figures.
