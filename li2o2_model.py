"""

Li-O2 Battery Model:
    This model examines the reactions taking place within the carbon-based
    cathode of a Li-O2 battery. Electrolyte = 1 M LiTFSI in TEGDME

"""

""" Load any needed modules """
"============================================================================"
# Brings in other complex commands as shortcuts
# Shortcut so that code doesn't need to be written out again
import numpy as np      # Support for multidimensional arrays and functions
import cantera as ct    # Open source chemistry toolbox
import matplotlib.pyplot as plt    # Plotting functions
from scipy.integrate import solve_ivp    #Integrator

""" Read user inputs and initialize variables, vectors, etc. """
"============================================================================"
from li2o2_init import objs, params, SVptr, pltptr, SV_0, tspan

""" Import derivative/residual function """
"""  TODO: THIS NEEDS TO BE MOVED OUT INTO ITS OWN FUNCTION, SIMILAR TO 'INIT' ABOVE"""
# Define function to solve
def LiO2_func(t,SV,params,objs,SVptr):
    # initializes derivative to all 0, allocates memory:
    dSVdt = np.zeros_like(SV)

    dPhidt = np.zeros_like(SV)
    dRhoOxidedt = np.zeros_like(SV)
    dRhoElytedt = np.zeros_like(SV)

    #  Pulls phase objects out of storage in 'objs' so they can be used more #    conveniently.
    gas = objs['gas']
    ca_bulk = objs['ca_bulk']
    elyte = objs['elyte']
    oxide = objs['oxide']
    ca_surf = objs['ca_surf']
    air_elyte = objs['air_elyte']
    Li_bulk = objs['Li_bulk']
    Li_surf = objs['Li_surf']

    # Initialize electronic and ionic currents
    i_ext = params['i_ext']     # [A]
    i_io = np.zeros(params['Ny'] + 1)     # initialize ionic current vector
    i_el = np.zeros(params['Ny'] + 1)     # initialize electronic current vector
    i_el[0] = i_ext             # electric current at air/cathode boundary
    i_io[-1] = i_ext            # ionic current at cathode/elyte boundary

    # Initialize electrolyte species flux vectors
    J_k_elyte = np.zeros((params['Ny']+1, int(elyte.n_species)))
    J_k_elyte[-1,params['i_Li_elyte']] = i_ext/ct.faraday

    # Set electric potentials
    # Are these constant values of 0?  SCD - Yes.  We assume that the cathode 
    #   is at an electric potential of zero (i.e. it is the reference potential)
    Phi_elyte = SV[SVptr['phi_dl']]
    ca_bulk.electric_potential = 0
    oxide.electric_potential = 0
    elyte.electric_potential = Phi_elyte

    # Look at first node:
    j = 0

    # Oxide volume fraction:
    eps_oxide = SV[SVptr['rho oxide'][j]] / oxide.density_mass    
    # Electrolyte volume fraction:
    eps_elyte = params['eps_elyte_0'] - (eps_oxide - params['eps_oxide_0'])
    # Electrolytte mass density (kg per m3 of electrolyte phase)
    rho_elyte = (sum(SV[SVptr['rho_k elyte'][j]])) / eps_elyte
    # Set electrolyte state.  Species mass densities are normalized by Cantera:
    elyte.TP = params['T'], ct.one_atm
    elyte.Y = SV[SVptr['rho_k elyte'][j]]

    # Calculate net production rates at interface
    sdot = ca_surf.net_production_rates                 # interface production rates

    # Calculate Faradaic current
    i_far = -ct.faraday*ca_surf.get_net_production_rates(ca_bulk)            # Faradaic current

    # Calculate change in oxide concentration
    W_oxide = oxide.mean_molecular_weight             # oxide molecular weight
    A_int_avail = params['A_int'] - eps_oxide / params['th_oxide']  # available interface area on carbon particle

    sdot_oxide = ca_surf.get_net_production_rates(oxide)
    dRhoOxidedt = sdot_oxide * A_int_avail * W_oxide

    # Calculate change in double layer potential
    i_dl = (i_io[0] - i_io[-1])*params['dyInv'] - i_far*A_int_avail    # double layer current
    dPhidt = i_dl / (params['C_dl']*params['A_int'])                           # double layer potential

    # Calculate change in electrolyte concentrations
    J_k_elyte_in = J_k_elyte[j,:]
    J_k_elyte_out  = J_k_elyte[j+1,:]
    sdot_elyte = ca_surf.get_net_production_rates(elyte)

    
    dRhoElytedt = (J_k_elyte_in - J_k_elyte_out)*params['dyInv'] + \
        (sdot_elyte * A_int_avail * elyte.molecular_weights)

    # Load differentials into dSVdt
    # Change in time for below variables
    dSVdt[SVptr['phi_dl']] = dPhidt                            # double layer potential
    dSVdt[SVptr['rho oxide']] = dRhoOxidedt                     # oxide concentration
    dSVdt[SVptr['rho_k elyte']] = dRhoElytedt                     # electrolyte concentration

    if 0:
        print('-----------------------------------------------------------')
        print('Phi =',Phi_cathode)
        print('eps_oxide =',eps_oxide)
        print('eps_elyte =',eps_elyte)
        print('sdot =',sdot)
        print('i_Far =',i_far)
        print('i_dl =',i_dl)
        elyte()

    return dSVdt

# Solve function using IVP solver
# Solves variables with different changes in time?
SV = solve_ivp(lambda t, y: LiO2_func(t,y,params,objs,SVptr), [0, tspan], SV_0, method='BDF',atol=params['atol'],rtol=params['rtol'])

#    Phi_dl = SV.y[SVptr['phi_dl'],-1]

#    return SV

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
