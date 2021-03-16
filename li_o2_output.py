""" li_o2_output.py

    Prepare and save figures, organize and save data
"""
""" Plot solutions to concentrations and potentials """

from matplotlib import pyplot as plt
"============================================================================"
def plot_profiles(SV, SVptr, objs, params, pltptr):
    legends = []
    [legends.append(str(i+1)) for i in range(params['N_y'])]

    fig, axs = plt.subplots(4,2, sharex=True, figsize=(10., 12.))
    fig.tight_layout

    # Plot Cell Potential:
    [axs[0,0].plot(SV.t, -SV.y[SVptr['phi_dl']][j]) for j in 
        range(params['N_y'])]
    axs[0,0].set(xlabel='Time (s)', ylabel='Double Layer Potential (V)')
    fig.legend(legends)

    # Plot Li2O2 volume fraction:
    [axs[1,0].plot(SV.t, SV.y[SVptr['eps oxide'][j]]) for j in 
        range(params['N_y'])]
    axs[1,0].set(xlabel='Time (s)', ylabel='Oxide volume fraction')

    oxide = objs['oxide']
    elyte = objs['elyte']
    eps_oxide = SV.y[SVptr['eps oxide'][:]]    # oxide volume fraction
    eps_elyte = params['eps_elyte_0'] - (eps_oxide - params['eps_oxide_0'])
    A_int_avail = params['A_int'] - eps_oxide / params['th_oxide']

    [axs[2,0].plot(SV.t,eps_elyte[j]) for j in range(params['N_y'])]
    axs[2,0].set(xlabel='Time (s)', ylabel='Elyte Volume Fraction')

    [axs[1,1].plot(SV.t, SV.y[SVptr['rho_k elyte'][j,2]]/eps_elyte[j]) for j in 
        range(params['N_y'])]
    axs[1,1].set(xlabel='Time (s)', ylabel=elyte.species_names[2]+
        ' kg/m3 elyte')

    [axs[2,1].plot(SV.t, SV.y[SVptr['rho_k elyte'][j,4]]/eps_elyte[j]) for j in 
        range(params ['N_y'])]
    axs[2,1].set(xlabel='Time (s)', ylabel=elyte.species_names[4]+
        ' kg/m3 elyte')

    [axs[3,0].plot(SV.t,A_int_avail[j,:]) for j in range(params['N_y'])]
    axs[3,0].set(xlabel='Time (s)', ylabel='Available Area (m2)')

    [axs[0,1].plot(SV.t, SV.y[SVptr['rho_k elyte'][j,2]]) for j in range(params
        ['N_y'])]
    axs[0,1].set(xlabel='Time (s)', ylabel=elyte.species_names[2]+
        ' kg/m3 total')

    axs[3,1].plot(SV.t, SV.y[pltptr['O2']]/eps_elyte[0], SV.t, 
        SV.y[pltptr['Li+']]/eps_elyte[0], SV.t, 
        SV.y[pltptr['PF6-']]/eps_elyte[0], SV.t, 
        SV.y[pltptr['EC']]/eps_elyte[0], SV.t, SV.y[pltptr['EMC']]/eps_elyte[0])
    axs[3,1].legend(['O2','Li+','PF6-','EC','EMC'], loc="upper right")
    axs[3,1].set(xlabel='Time (s)', ylabel='Concentration (kg/m3 tot)')

    fig.tight_layout()
    plt.show()