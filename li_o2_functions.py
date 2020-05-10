import numpy as np      # Support for multidimensional arrays and functions
import cantera as ct    # Open source chemistry toolbox

def read_state(SV, SVptr, j):
    # double layer (i.e. elyte) electric potential:
    phi_elyte = SV[SVptr['phi_dl'][j]]

    # Oxide volume fraction:
    eps_oxide = SV[SVptr['eps oxide'][j]]

    # Electrolyte species mass fracionts:
    rho_k_elyte = SV[SVptr['rho_k elyte'][j]]

    return phi_elyte, eps_oxide, rho_k_elyte


# Define function for the integrator/solver:
def li_o2_residual(t,SV,params,objs,SVptr):

    #  Pulls phase objects out of storage in 'objs' so they can be used more #    conveniently.
    gas = objs['gas']
    ca_bulk = objs['ca_bulk']
    elyte = objs['elyte']
    oxide = objs['oxide']
    ca_surf = objs['ca_surf']
    air_elyte = objs['air_elyte']
    Li_bulk = objs['Li_bulk']
    Li_surf = objs['Li_surf']

    # Faraday's constant (C/kmol e-)
    F = ct.faraday
    # oxide molar volume
    Vbar_oxide = oxide.mean_molecular_weight / oxide.density_mass       

    R_io = 1e6

    # initializes derivative to all 0, allocates memory:
    dSVdt = np.zeros_like(SV)

    # Initialize electronic and ionic currents
    i_ext = params['i_ext']     # [A]
    i_io = np.zeros(params['N_y'] + 1)     # initialize ionic current vector
    i_el = np.zeros(params['N_y'] + 1)     # initialize electronic current vector
    i_el[0] = i_ext             # electric current at air/cathode boundary
    i_io[-1] = i_ext            # ionic current at cathode/elyte boundary

    # Initialize electrolyte species flux vectors
    J_k_elyte = np.zeros((params['N_y']+1, int(elyte.n_species)))
    J_k_elyte[-1,params['i_Li_elyte']] = i_ext/F

    j = 0
    phi_elyte, eps_oxide, rho_k_elyte = read_state(SV, SVptr, j)
    # Set electric potentials
    # Are these constant values of 0?  SCD - Yes.  We assume that the cathode 
    #   is at an electric potential of zero (i.e. it is the reference potential)
    ca_bulk.electric_potential = 0
    oxide.electric_potential = 0
    elyte.electric_potential = phi_elyte 

    oxide.TP = params['T'], ct.one_atm
    elyte.TP = params['T'], ct.one_atm
    #Species mass densities are normalized by Cantera to give mass fractions:
    elyte.Y = rho_k_elyte

    for j in range(params['N_y']-1):
        
        # Calculate net production rates at interface
        #    Oxide species:
        sdot_oxide = ca_surf.get_net_production_rates(oxide)
        #    Electrolyte species:
        sdot_elyte = ca_surf.get_net_production_rates(elyte)
        
        # Calculate Faradaic current
        i_far = -F*ca_surf.get_net_production_rates(ca_bulk) # Faradaic current

        # Calculate change in oxide concentration
        A_int_avail = params['A_int'] - eps_oxide / params['th_oxide']  # available interface area on carbon particle

        # Read out state from next node and set Cantera object states:
        phi_elyte_next, eps_oxide_next, rho_k_elyte_next = read_state(SV, SVptr, j+1)

        # SCDm 09 May 2020--This is backwards, for now, just to implement a 
        #  simple transport function.  The eventual implementation should 
        #  calculate the species fluxes, and use this to calculate i_io.
        i_io[j+1] = (phi_elyte - phi_elyte_next)*params['dyInv']/R_io
        J_k_elyte[j,params['i_Li_elyte']] = i_io[j]/F

        dEpsOxide_dt = sdot_oxide * A_int_avail * Vbar_oxide   

        # Double layer current
        i_dl = (i_io[j] - i_io[j+1])*params['dyInv'] - i_far*A_int_avail 
        # Change in double layer potential per unit time 
        dPhi_dt = i_dl / (params['C_dl']*params['A_int'])

        # Calculate change in electrolyte concentrations
        J_k_elyte_in = J_k_elyte[j,:]
        J_k_elyte_out  = J_k_elyte[j+1,:]
        
        dRhoElyte_dt = (J_k_elyte_in - J_k_elyte_out)*params['dyInv'] + \
            (sdot_elyte * A_int_avail * elyte.molecular_weights)

        # Load differentials into dSVdt
        # Change in time for below variables
        dSVdt[SVptr['phi_dl'][j]] = dPhi_dt                            # double layer potential
        dSVdt[SVptr['eps oxide'][j]] = dEpsOxide_dt                     # oxide concentration
        dSVdt[SVptr['rho_k elyte'][j]] = dRhoElyte_dt                     # electrolyte concentration

        phi_elyte = phi_elyte_next
        eps_oxide = eps_oxide_next
        rho_k_elyte = rho_k_elyte_next
        # Set electric potentials
        ca_bulk.electric_potential = 0
        oxide.electric_potential = 0
        elyte.electric_potential = phi_elyte 

        oxide.TP = params['T'], ct.one_atm
        elyte.TP = params['T'], ct.one_atm
        #Species mass densities are normalized by Cantera to give mass fractions:
        elyte.Y = rho_k_elyte
        
        
    j = params['N_y']-1
    print(i_io)  
    # Calculate net production rates at interface
    #    Oxide species:
    sdot_oxide = ca_surf.get_net_production_rates(oxide)
    #    Electrolyte species:
    sdot_elyte = ca_surf.get_net_production_rates(elyte)

    # Calculate Faradaic current
    i_far = -F*ca_surf.get_net_production_rates(ca_bulk) # Faradaic current

    # Calculate change in oxide concentration
    A_int_avail = params['A_int'] - eps_oxide / params['th_oxide']  # available interface area on carbon particle

    dEpsOxide_dt = sdot_oxide * A_int_avail * Vbar_oxide  

    # Calculate change in double layer potential
    i_dl = (i_io[j] - i_io[j+1])*params['dyInv'] - i_far*A_int_avail    # double layer current
    dPhi_dt = i_dl / (params['C_dl']*params['A_int'])                           # double layer potential

    # Calculate change in electrolyte concentrations
    J_k_elyte_in = J_k_elyte[j,:]
    J_k_elyte_out  = J_k_elyte[j+1,:]
    
    dRhoElyte_dt = (J_k_elyte_in - J_k_elyte_out)*params['dyInv'] + \
        (sdot_elyte * A_int_avail * elyte.molecular_weights)

    # Load differentials into dSVdt
    # Change in time for below variables
    dSVdt[SVptr['phi_dl'][j]] = dPhi_dt                            # double layer potential
    dSVdt[SVptr['eps oxide'][j]] = dEpsOxide_dt                     # oxide concentration
    dSVdt[SVptr['rho_k elyte'][j]] = dRhoElyte_dt                     # electrolyte concentration

    return dSVdt