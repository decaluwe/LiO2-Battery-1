import numpy as np      # Support for multidimensional arrays and functions
import cantera as ct    # Open source chemistry toolbox

# Define function for the integrator/solver:
def li_o2_residual(t,SV,params,objs,SVptr):
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
    i_io = np.zeros(params['N_y'] + 1)     # initialize ionic current vector
    i_el = np.zeros(params['N_y'] + 1)     # initialize electronic current vector
    i_el[0] = i_ext             # electric current at air/cathode boundary
    i_io[-1] = i_ext            # ionic current at cathode/elyte boundary

    # Initialize electrolyte species flux vectors
    J_k_elyte = np.zeros((params['N_y']+1, int(elyte.n_species)))
    J_k_elyte[-1,params['i_Li_elyte']] = i_ext/ct.faraday

    # Set electric potentials
    # Are these constant values of 0?  SCD - Yes.  We assume that the cathode 
    #   is at an electric potential of zero (i.e. it is the reference potential)
    Phi_elyte = SV[SVptr['phi_dl']]
    ca_bulk.electric_potential = 0
    oxide.electric_potential = 0
    elyte.electric_potential = Phi_elyte

    for j in range(params['N_y']):
        
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
        dSVdt[SVptr['phi_dl'][j]] = dPhidt                            # double layer potential
        dSVdt[SVptr['rho oxide'][j]] = dRhoOxidedt                     # oxide concentration
        dSVdt[SVptr['rho_k elyte'][j]] = dRhoElytedt                     # electrolyte concentration

    return dSVdt