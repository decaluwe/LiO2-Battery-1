"""

Li-O2 Battery Model:
    This model examines the reactions taking place within the carbon-based
    cathode of a Li-O2 battery. Electrolyte = 1 M LiTFSI in TEGDME

"""

""" Load any needed modules """
"============================================================================"
from scipy.integrate import solve_ivp   # Integrator

""" Read user inputs and initialize variables, vectors, etc. """
"============================================================================"
# Import cantera objects, parameters, pointers, initial solution vector SV_0, 
#    and residual function
from li_o2_init import objs, params, SVptr, pltptr, SV_0, tspan, li_o2_residual
from li_o2_terminate import voltage_min

flag_discharge, flag_charge = False, False

# Possible steps include 'Equilibrating', 'Discharging', and 'Charging.'
steps = params['n_cycles']*('Equilibrating', 'Discharging')
currents = params['n_cycles']*([-1e-16, params['i_ext'], 2e-12, 
    -params['i_ext']])

# Print a blank line:
print('\n')

for i, step in enumerate(steps):
    print(step,'...\n')
    params['i_ext'] = currents[i]
    print('     Current = ', round(currents[i],3),'\n')
    if step=='Discharging':
        flag_discharge = True
        voltage_min.terminal = True
        SV_discharge =  solve_ivp(li_o2_residual, [0, tspan], SV_0, 
            method='BDF', args=(params,objs,SVptr), events=voltage_min, atol=params['atol'],rtol=params['rtol'])
        SV_0 = SV_discharge.y[:,-1]
    elif step=='Charging':
        flag_charge = True
        voltage_min.terminal = True
        SV_charge =  solve_ivp(li_o2_residual, [0, tspan], SV_0, method='BDF',
            args=(params, objs, SVptr), atol=params['atol'],rtol=params['rtol'])
        SV_0 = SV_charge.y[:,-1]
    else:
        voltage_min.terminal = False
        SV_equil = solve_ivp(li_o2_residual, [0, tspan], SV_0, method='BDF', 
            args=(params, objs, SVptr), events=voltage_min, atol=params['atol'],rtol=params['rtol'])
        SV_0 = SV_equil.y[:,-1]

print('Done with simulation. Preparing outputs.\n')
from li_o2_output import plot_profiles
import matplotlib.pyplot as plt

# Plot discharge profiles:
if flag_discharge:
    plot_profiles(SV_discharge, SVptr, objs, params, pltptr)
# plot_profiles(SV_equil, SVptr, objs, params, pltptr)
# Plot charge profiles:
if flag_charge:
    plot_profiles(SV_charge, SVptr, objs, params, pltptr)
plt.show()