"========== USER INPUTS for 1D Li-O2 Battery simulation =========="

phi_elyte_init = -3.19            # double layer voltage [V]
eps_elyte_init = 0.5                # initial electrolyte volume fraction [-]
eps_oxide_init = 1e-12              # initial oxide volume fraction [-]
eps_binder_init = 0.                # initial binder volume fraction [-]
eps_carbon = 1. - eps_elyte_init - eps_binder_init - eps_oxide_init      # initial carbon volume fraction [-]

# Tolerances(absolute and relative)
atol = 1e-10
rtol = 2.5e-6

i_ext = -1e1              # [A/m2]    c-rate will calculate this
cap = 1e-3*2.1733333       # battery capacity

# How long integration goes [s]. Should eventually pick a c-rate and sub in
tspan = 230450000*abs(1e-3/i_ext)#3100604#4305 #7824

N_x = 1                      # 1D model track state in each cell and each layer
N_y = 5                      # no. of cells in the y-direction
Nvars = 3                   # no. of variables
th_ca = 50e-6               # cathode thickness [m]
d_part = 10e-6              # carbon particle diameter [m]
d_oxide = 2e-6              # oxide particle diameter [m]
th_oxide = 5e-6             # thickness of oxide ellipsoid [m]
C_dl = 1.1e-6               # double layer capacitance [F/m2]

TP = 300, 101325             # inital temp, pressure [K, Pa]

# Electrolyte bulk diffusion coefficients
D_k_elyte = {}
D_k_elyte['Li+(e)'] = 4e-11          # bulk diff coeff Li+ in elyte (m2/s)
D_k_elyte['PF6-(e)'] = 4e-13         # bulk diff coeff PF6- in elyte (m2/s)
D_k_elyte['O2(e)'] = 7e-12           # bulk diff coeff O2 in elyte (m2/s)
D_k_elyte['C3H4O3(e)'] = 1           # EC diffusion is fast
D_k_elyte['C4H8O3(e)'] = 1           # EMC diffusion is fast

# Bruggeman coeff. for tortuosity factor: tau_fac = eps_elyte^(-n_bruggeman)
n_bruggeman = 0.5 

ctifile = 'li_o2_battery.yaml'     # Cantera input file
li_elyte_name = 'Li+(e)'

" Specifyt the residual function.  Import it as 'li_o2_residual"
" This function should be written in li_o2_functions--the only term"
" to edit is therefore the middle name."
"============================================================================"
from li_o2_functions import li_o2_residual as li_o2_residual