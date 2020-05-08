"""

Li-O2 Battery Model:
    This model examines the reactions taking place within the carbon-based
    cathode of a Li-O2 battery. Electrolyte = 1 M LiTFSI in TEGDME

"""

""" Load any needed modules """
"============================================================================"
# Brings in other complex commands as shortcuts
# Shortcut so that code doesn't need to be written out again
import matplotlib.pyplot as plt    # Plotting functions
from scipy.integrate import solve_ivp    #Integrator

""" Read user inputs and initialize variables, vectors, etc. """
"============================================================================"
from li_o2_init import objs, params, SVptr, pltptr, SV_0, tspan, li_o2_residual

# Solve function using IVP solver
SV = solve_ivp(lambda t, y: li_o2_residual(t,y,params,objs,SVptr), [0, tspan], \
     SV_0, method='BDF',atol=params['atol'],rtol=params['rtol'])


""" Plot solutions to concentrations and potentials """
"============================================================================"
plt.figure(1)
plt.plot(SV.t,-SV.y[SVptr['phi_dl'][0]])
plt.xlabel('Time (s)')
plt.ylabel('Double Layer Potential (V)')

#plt.figure(2)
#plt.plot(SV.t,SV.y[SVptr['oxide']])
#plt.xlabel('Time (s)')
#plt.ylabel('Oxide Concentration (kg/m3)')

oxide = objs['oxide']

eps_oxide = SV.y[SVptr['rho oxide']] / oxide.density_mass      # oxide volume fraction
eps_elyte = params['eps_elyte_0'] - (eps_oxide - params['eps_oxide_0'])
A_int_avail = params['A_int'] - eps_oxide / params['th_oxide']

#plt.figure(3)
#plt.plot(SV.t,eps_elyte)
#plt.xlabel('Time (s)')
#plt.ylabel('Elyte Volume Fraction')
#plt.show()

#plt.figure(4)
#plt.plot(SV.t/3600 * -i_ext,SV.y[SVptr['phi_dl']])
#plt.xlabel('Capacity (Ah/m2)')
#plt.ylabel('Voltage (V)')

plt.figure(5)
plt.plot(SV.t,A_int_avail[0,:])
plt.xlabel('Time (s)')
plt.ylabel('Available Area (m2)')

plt.show()

#plt.figure(3)
#plt.plot(SV.t,SV.y[pltptr['O2']],SV.t,SV.y[pltptr['Li+']],SV.t,SV.y[pltptr['PF6-']],SV.t,SV.y[pltptr['EC']],SV.t,SV.y[pltptr['EMC']])
#plt.legend(['O2','Li+','PF6-','EC','EMC'])
#plt.xlabel('Time (s)')
#plt.ylabel('Electrolyte Concentration (kg/m3)')
#plt.show()


#t = SV.t
#dl = SV.y[SVptr['phi_dl']]
#Ck_ox = SV.y[SVptr['oxide']]
#
#df = DataFrame({'Time': t, 'Double Layer': dl, 'Oxide Concentration': Ck_ox})
#
#with ExcelWriter('path_to_file.xlsx') as writer:
#    df.to_excel(writer)
