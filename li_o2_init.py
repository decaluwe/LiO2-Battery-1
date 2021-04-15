" Initialize variables, parameters, and data structures for simulation"

" Load necessary modules "
"============================================================================"
import numpy as np
import cantera as ct

" Load in user inputs "
"============================================================================"
from li_o2_inputs import *

" Use inputs to calculate relevant quantities:"
"============================================================================"
# Geometric calculations
V_part = 4./3. * np.pi * (d_part / 2)**3  # particle volume [m3]
A_part = 4. * np.pi * (d_part / 2)**2    # particle surface area [m2]
A_int = eps_carbon * A_part / V_part      # interface area [m2/m3 total]
A_oxide = np.pi * d_oxide**2 / 4.        # oxide area contacting carbon particle
V_oxide = 2./3. * np.pi * (d_oxide/2)**2 * th_oxide   # oxide volume [m3]

# Import necessary phases from Cantera
# Make objects to handle calculations
gas = ct.Solution(ctifile,'oxygen')
ca_bulk = ct.Solution(ctifile,'conductor')
elyte = ct.Solution(ctifile,'electrolyte')
oxide = ct.Solution(ctifile,'li-oxide')
ca_surf = ct.Interface(ctifile,'cathode-surf',[elyte, oxide, ca_bulk])
air_elyte = ct.Interface(ctifile,'air-elyte-interface',[gas, elyte])
Li_conductor  = ct.Solution(ctifile,'conductor')
Li_bulk = ct.Solution(ctifile,'lithium-bulk')
Li_surf = ct.Interface(ctifile,'anode-surf',[Li_bulk, elyte, Li_conductor])

# Sets states for the objects
oxide.TP = TP
elyte.TP = TP
ca_surf.TP = TP
ca_bulk.TP = TP

# Store these phases in a common 'objs' dict.
# This takes the data from Cantera and stores data for use in a Python 
# dictionary so we can efficiently pass them to various functions. A dictionary
# associates a key word with a value.
objs = {}
objs['gas'] = gas
objs['ca_bulk'] = ca_bulk
objs['elyte'] = elyte
objs['oxide'] = oxide
objs['ca_surf'] = ca_surf
objs['air_elyte'] = air_elyte
objs['Li_bulk'] = Li_bulk
objs['Li_surf'] = Li_surf

# Battery capacity: Ah/m2
Capacity = th_ca*(1-eps_carbon)/oxide.partial_molar_volumes[0]*2*ct.faraday/3600

if C_rate:
    i_ext = -C_rate*Capacity
else:
    C_rate = i_ext/Capacity

tspan = 3600/C_rate
# Store parameters in a common 'params' dict
# This allows for an unknown amount of arguments and any type of data to be stored at once?
# Dictionary to keep like things together
params = {}
params['i_ext'] = i_ext
params['n_cycles'] = n_cycles
params['T'] = TP[0]
params['eps_elyte_0'] = eps_elyte_init
params['eps_oxide_0'] = eps_oxide_init
params['rtol'] = rtol
params['atol'] = atol
params['N_y'] = int(N_y)
params['A_int'] = A_int
params['th_oxide'] = th_oxide
params['dyInv'] = N_y/th_ca
params["C_dl"] = C_dl
params['i_Li_elyte'] = elyte.species_index(li_elyte_name)
params['save_name'] = save_name

# Bulk electrolyte diffusion coefficients:
params['D_o_elyte'] = np.zeros_like(elyte.X)
params['z_k_elyte'] = np.zeros_like(elyte.X)
for j in D_k_elyte:
    params['D_o_elyte'][elyte.species_index(j)] = D_k_elyte[j]
    params['z_k_elyte'][elyte.species_index(j)] = elyte.species(j).charge
    
params['n bruggeman'] = n_bruggeman
# Store solution vector pointers in a common 'SVptr' dict
# One spot for each of the species, can change Cantera file to change # and automatically adapts
SVptr = {}
# Variables per finite volume (aka 'node'): elyte electric potential, oxide 
#     density, electrolyte species densities
nvars_node = int(elyte.n_species + 2)

# double layer potential in solution vector SV
SVptr['phi_dl'] = np.arange(0,nvars_node*N_y,nvars_node, \
    dtype='int')                                     
# oxide volume fraction in solution vector SV
SVptr['eps oxide'] = np.arange(1,nvars_node*N_y, nvars_node, \
    dtype='int')                                      
# electrolyte species mass densities in solution vector SV
SVptr['rho_k elyte'] = np.ndarray(shape=(N_y, elyte.n_species),\
    dtype='int')       
for i in range(N_y):
    SVptr['rho_k elyte'][i,:] = range(2 + i*nvars_node, 2 + i*nvars_node + \
         elyte.n_species)




# Store plot pointers in a common 'pltptr' dict
# Stored values to be used in plot
pltptr = {}
pltptr['O2'] = 6
pltptr['Li+'] = 4
pltptr['PF6-'] = 5
pltptr['EC'] = 2
pltptr['EMC'] = 3

# Set inital values
# This shows the equations for change over time?
eps_oxide_init = params['eps_oxide_0']          # oxide concentraion
rho_elyte_init = elyte.Y*elyte.density*params['eps_elyte_0']  # electrolyte concentrations
SV0 = np.r_[phi_elyte_init, eps_oxide_init, rho_elyte_init]   # store in an array
SV_0 = np.tile(SV0, N_y)                   # tile SV0 based on discretization