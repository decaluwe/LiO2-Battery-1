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

# Save user-defined i_ext:
i_input = params['i_ext']

steps = ('Equilibrating', 'Discharging')#, 'Equilibrating', 'Charging')
currents = ([-1e-16, params['i_ext'], 1e-16, -params['i_ext']])

# Clear a blank line:
print('\n')

for i, step in enumerate(steps):
    print(step,'...\n')
    params['i_ext'] = currents[i]
    print('     Current = ', round(currents[i],3),'\n')
    if step=='Discharging':
        SV_discharge =  solve_ivp(lambda t, y: li_o2_residual(t,y,params,objs,  
            SVptr), [0, tspan], SV_0, method='BDF',atol=params['atol'],rtol=params['rtol'])
        SV_0 = SV_discharge.y[:,-1]
    elif step=='Charging':
        SV_charge =  solve_ivp(lambda t, y: li_o2_residual(t,y,params,objs,  
            SVptr), [0, tspan], SV_0, method='BDF',atol=params['atol'],rtol=params['rtol'])
        SV_0 = SV_charge.y[:,-1]
    else:
        SV_equil = solve_ivp(lambda t, y: li_o2_residual(t,y,params,objs,SVptr),
            [0, tspan], SV_0, method='BDF',atol=params['atol'],rtol=params['rtol'])
        SV_0 = SV_equil.y[:,-1]

print('Done with simulation. Preparing outputs.\n')
from li_o2_output import plot_profiles

# Plot discharge profiles:
plot_profiles(SV_discharge, SVptr, objs, params, pltptr)