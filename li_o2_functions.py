import numpy as np      # Support for multidimensional arrays and functions
import cantera as ct    # Open source chemistry toolbox

# Define function for the integrator/solver:
def li_o2_residual(t,SV,params,objs,SVptr):

    #  Pulls phase objects out of storage in 'objs' so they can be used more #    conveniently.
    gas, ca_bulk, elyte, oxide, ca_surf, air_elyte, Li_bulk, Li_surf = \
        read_cantera_objs(objs)
    
    # Read out and store some constants:
    # Faraday's constant (C/kmol e-)
    F = ct.faraday
    # oxide molar volume
    Vbar_oxide = oxide.mean_molecular_weight / oxide.density_mass       
    # ELectrolyte molar weights
    W_elyte = elyte.molecular_weights

    # TEMPORARY:
    R_io = 3e6

    # Initializes derivative to all 0, allocates memory:
    dSVdt = np.zeros_like(SV)

    # Initialize electronic and ionic currents
    i_ext = params['i_ext']     # [A]
    i_io = np.zeros(params['N_y'] + 1)     # ionic current vector
    i_el = np.zeros(params['N_y'] + 1)     # electronic current vector
    i_el[0] = i_ext             # electric current at air/cathode boundary
    i_io[-1] = i_ext            # ionic current at cathode/elyte boundary

    # Initialize electrolyte species flux vectors
    J_k_elyte = np.zeros((params['N_y']+1, int(elyte.n_species)))
    J_k_elyte[-1,params['i_Li_elyte']] = W_elyte[params['i_Li_elyte']]*i_ext/F

    # Electrolyte species production rates due to double layer current:
    sdot_dl = np.zeros_like(elyte.X)

    # Look at node adjacent to air/cathode interface:
    j = 0
    phi_elyte, eps_oxide, rho_k_elyte = read_state(SV, SVptr, j)
    eps_elyte = params['eps_elyte_0'] - (eps_oxide - params['eps_oxide_0'])
    # Set cantera object states:
    ca_bulk, oxide, elyte, air_elyte, gas = set_state(phi_elyte, rho_k_elyte, \
        params, ca_bulk, oxide, elyte, air_elyte, gas) 

    # Set the flux into the electrolyte at the air/electrolyte interface equal 
    #    to the reaction rate:
    J_k_elyte[0,:] = air_elyte.get_net_production_rates(elyte)*W_elyte

    for j in range(params['N_y']-1):
        
        # Calculate species net production rates from Cantera objects:
        sdot_oxide, sdot_elyte, i_far = read_rates(oxide, elyte, ca_bulk, \
            ca_surf)
        
        # Calculate change in oxide concentration
        A_int_avail = params['A_int'] - eps_oxide / params['th_oxide']  # available interface area on carbon particle

        # Read out state from next node and set Cantera object states:
        phi_elyte_next, eps_oxide_next, rho_k_elyte_next = read_state(SV,\
            SVptr, j+1)
        eps_elyte_next = params['eps_elyte_0'] - (eps_oxide_next -\
            params['eps_oxide_0'])
        
        J_k_elyte[j+1, :], i_io_check = mass_fluxes(params, W_elyte, \
            rho_k_elyte, rho_k_elyte_next, phi_elyte, phi_elyte_next, \
            eps_elyte,  eps_elyte_next)
        
        # print(i_io_check)
        # SCDm 09 May 2020--This is backwards, for now, just to implement a 
        #  simple transport function.  The eventual implementation should 
        #  calculate the species fluxes, and use this to calculate i_io.
        i_io[j+1] = i_io_check#(phi_elyte - phi_elyte_next)*params['dyInv']/R_io
        # J_k_elyte[j+1,params['i_Li_elyte']] = i_io[j+1]/F

        dEpsOxide_dt = sdot_oxide * A_int_avail * Vbar_oxide   

        # Double layer current
        i_dl = (i_io[j] - i_io[j+1])*params['dyInv'] - i_far*A_int_avail 
        # Change in double layer potential per unit time 
        dPhi_dt = i_dl / (params['C_dl']*params['A_int'])
        # print(j, i_far, sdot_elyte[2], i_dl)
        # Double layer current consumes Li+(e)
        sdot_dl[params['i_Li_elyte']] = -i_dl/F

        # Calculate change in electrolyte concentrations       
        dRhoElyte_dt = ((J_k_elyte[j,:] - J_k_elyte[j+1,:])*params['dyInv'] + \
            sdot_elyte*A_int_avail*W_elyte + sdot_dl*W_elyte) - \
            rho_k_elyte*dEpsOxide_dt
        
        # Load differentials into dSVdt
        # Change in time for below variables
        dSVdt[SVptr['phi_dl'][j]] = dPhi_dt                            # double layer potential
        dSVdt[SVptr['eps oxide'][j]] = dEpsOxide_dt                     # oxide concentration
        dSVdt[SVptr['rho_k elyte'][j]] = dRhoElyte_dt                     # electrolyte concentration

        # Prepare for next 'j' value: re-assign "next" to "current" variables:
        phi_elyte = phi_elyte_next
        eps_oxide = eps_oxide_next
        rho_k_elyte = rho_k_elyte_next
        
        # Set cantera object states:
        ca_bulk, oxide, elyte, air_elyte, gas = set_state(phi_elyte, \
            rho_k_elyte, params, ca_bulk, oxide, elyte, air_elyte, gas)        
        
    
    # Calculations for the node adjacent to the electrolyte separator:
    j = params['N_y']-1
     
    # Calculate species net production rates from Cantera objects:
    sdot_oxide, sdot_elyte, i_far = read_rates(oxide, elyte, ca_bulk, ca_surf)

    # Available interface area on carbon particle
    A_int_avail = params['A_int'] - eps_oxide / params['th_oxide']  
    
    # Calculate change in oxide concentration
    dEpsOxide_dt = sdot_oxide * A_int_avail * Vbar_oxide  

    # Calculate change in double layer potential:
    #   Double layer current:
    i_dl = (i_io[j] - i_io[j+1])*params['dyInv'] - i_far*A_int_avail 
    
    # print(j, i_far, sdot_elyte[2], i_dl)
    # Double layer potential rate of change:
    dPhi_dt = i_dl / (params['C_dl']*params['A_int'])   
    # Convert the double layer current to an equivalent reaction rate consuming 
    #    Li+(e):   
    sdot_dl[params['i_Li_elyte']] = -i_dl/F

    # Calculate change in electrolyte concentrations   
    dRhoElyte_dt = ((J_k_elyte[j,:] - J_k_elyte[j+1,:])*params['dyInv'] + \
        sdot_elyte*A_int_avail*W_elyte + sdot_dl*W_elyte) - \
        rho_k_elyte*dEpsOxide_dt

    # print(j, J_k_elyte[j,2]*params['dyInv'], J_k_elyte[j+1,2]*params['dyInv'], \
    #     sdot_elyte[2]*A_int_avail*W_elyte[2], sdot_dl[2]*W_elyte[2])
    # Load differentials into dSVdt
    # Change in time for below variables
    dSVdt[SVptr['phi_dl'][j]] = dPhi_dt             # double layer potential
    dSVdt[SVptr['eps oxide'][j]] = dEpsOxide_dt     # oxide concentration
    dSVdt[SVptr['rho_k elyte'][j]] = dRhoElyte_dt   # electrolyte concentration
    
    return dSVdt

" ================= HELPER FUNCTIONS ================= "
    
def read_cantera_objs(objs):
    gas = objs['gas']
    ca_bulk = objs['ca_bulk']
    elyte = objs['elyte']
    oxide = objs['oxide']
    ca_surf = objs['ca_surf']
    air_elyte = objs['air_elyte']
    Li_bulk = objs['Li_bulk']
    Li_surf = objs['Li_surf']

    return gas, ca_bulk, elyte, oxide, ca_surf, air_elyte, Li_bulk, Li_surf

def read_state(SV, SVptr, j):
    # double layer (i.e. elyte) electric potential:
    phi_elyte = SV[SVptr['phi_dl'][j]]

    # Oxide volume fraction:
    eps_oxide = SV[SVptr['eps oxide'][j]]

    # Electrolyte species mass fracionts:
    rho_k_elyte = SV[SVptr['rho_k elyte'][j]]

    return phi_elyte, eps_oxide, rho_k_elyte

def set_state(phi_elyte, rho_k_elyte, params, ca_bulk, oxide, elyte, air_elyte,\
    gas):
    
    ca_bulk.electric_potential = 0
    oxide.electric_potential = 0
    elyte.electric_potential = phi_elyte 
    TP = params['T'], ct.one_atm

    oxide.TP = TP
    elyte.TP = TP
    gas.TP = TP
    air_elyte.TP = TP

    #Species mass densities are normalized by Cantera to give mass fractions:
    elyte.Y = rho_k_elyte

    return ca_bulk, oxide, elyte, air_elyte, gas

def read_rates(oxide, elyte, ca_bulk, ca_surf):
    sdot_oxide = ca_surf.get_net_production_rates(oxide)
    #    Electrolyte species:
    sdot_elyte = ca_surf.get_net_production_rates(elyte)
    
    # Calculate Faradaic current
    i_far = -ct.faraday*ca_surf.get_net_production_rates(ca_bulk) # Faradaic current

    return sdot_oxide, sdot_elyte, i_far 

def mass_fluxes(params, W_elyte, rho_k, rho_k_next, phi_elyte, phi_elyte_next, \
    eps_elyte, eps_elyte_next):
    
    # Store quantity F/RT:
    FoRT = ct.faraday/ct.gas_constant/params['T']
    
    # Current location:
    C_k = rho_k / W_elyte     # species molar densities
    Y_k = rho_k/sum(rho_k)    # species mass fractions
    X_k = C_k/sum(C_k)        # species mole fractions

    # Next location:
    C_k_next = rho_k_next / W_elyte        # species molar densities
    X_k_next = C_k_next/sum(C_k_next)      # species mole fractions

    # Take averages to find interface values.  Eventually this should be 
    #   weighted by the volume dimensions:
    C_k_elyte_int = 0.5*(C_k + C_k_next)
    eps_int = 0.5*(eps_elyte + eps_elyte_next)

    # Chemical diffusion and migration diffusion coefficients: 
    D_k_elyte = params['D_o_elyte']*eps_int**(1. + params['n bruggeman'])
    D_k_elyte_mig = D_k_elyte*params['z_k_elyte']*FoRT*C_k_elyte_int
    
    N_k = (D_k_elyte*(C_k/eps_elyte - C_k_next/eps_elyte_next)+ \
        D_k_elyte_mig*(phi_elyte - phi_elyte_next))*params['dyInv']

    i_io = np.dot(N_k, params['z_k_elyte'])*ct.faraday

    J_k = W_elyte*N_k

    return J_k, i_io