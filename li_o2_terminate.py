""" li_o2_terminate.py

    Determine events which would cause the simulation to terminate.
"""
import numpy as np

def voltage_min(t, SV, params, objs, SVptr):
    min_voltage = 5

    # Cathode double-layer voltage < -1.75 V
    for j in range(params['N_y']):
        min_voltage = min(min_voltage, -(SV[SVptr['phi_dl'][j]] + 2.0))

    return min_voltage