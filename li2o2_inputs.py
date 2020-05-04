"========== USER INPUTS for 1D Li-O2 Battery simulation =========="

phi_elyte_init = -3.19            # double layer voltage [V]
E_elyte_init = 0.5                # initial electrolyte volume fraction [-]
E_oxide_init = 1e-12              # initial oxide volume fraction [-]
E_binder_init = 0.                # initial binder volume fraction [-]
E_carbon = 1. - E_elyte_init - E_binder_init - E_oxide_init      # initial carbon volume fraction [-]

# Tolerances(absolute and relative)
atol = 1e-10
rtol = 2.5e-6

# How long integration goes [s]. Should eventually pick a c-rate and sub in
tspan = 7824

i_ext = -1e-3              # [A/m2]    c-rate will calculate this
cap = 1e-3*2.1733333       # battery capacity

Nx = 1                      # 1D model track state in each cell and each layer
Ny = 1                      # no. of cells in the y-direction
Nvars = 3                   # no. of variables
th_ca = 50e-6               # cathode thickness [m]
d_part = 10e-6              # carbon particle diameter [m]
d_oxide = 2e-6              # oxide particle diameter [m]
th_oxide = 5e-6             # thickness of oxide ellipsoid [m]
C_dl = 1.1e-6               # double layer capacitance [F/m2]

TP = 300, 101325             # inital temp, pressure [K, Pa]

ctifile = 'LiAir_mod.cti'     # Cantera input file