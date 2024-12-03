# This file calculates the efficiencies of a turbofan engine to then plot a sankey diagram
# Authors: Alfonso Rato, Jakub Bartmański 

def cycle_calculator(bypass, m_dot_c_isa, pi_fan, pi_hpc, pi_comb, T_04, pi_tot, eta_intake, eta_fan, eta_LPC, eta_HPC, eta_LPT, eta_HPT, eta_combustor, LHV):
    # Constants
    cp_air = 1000 
    kappa_air = 1.4
    cp_gas = 1150
    kappa_gas = 1.33

    M = 0.78
    T = 220
    P = 23842
    T_isa = 288.15
    P_isa = 101325

    # Calculate pressure ratio acroos LPC from total pressure ratio

    pi_lpc = pi_tot / (pi_fan * pi_hpc)

    # Calculate free stream velocity
    V = M * (kappa_air * R_air * T)**0.5

    # Calculate core massflow from corrected core massflow
    m_dot_c = m_dot_c_isa * P / P_isa * (T_isa / T)**0.5

    # Calculate bypass massflow from bypass ratio
    m_dot_f = bypass * m_dot_c

    # Calculate total temperature and pressure at free stream
    T_0_a = T * (1 + (kappa_air - 1) / 2 * M**2)
    P_0_a = P * (1 + (kappa_air - 1) / 2 * M**2)**(kappa_air / (kappa_air - 1))

    # Inlet 0 -> 2
    P_0_2 = P * (1 + eta_intake * (kappa_air - 1) / 2 * M**2)**(kappa_air / (kappa_air - 1))
    T_0_2 = T_0_a

    # Fan 2 -> 21
    P_0_21 = P_0_2 * pi_fan
    T_0_21 = T_0_2 * (1 + 1/eta_fan * (pi_fan**((kappa_air - 1) / kappa_air) - 1))

    # Bypass 2 -> 13
    P_0_13 = P_0_21
    T_0_13 = T_0_21

    # LPC 21 -> 25
    P_0_25 = P_0_21 * pi_lpc
    T_0_25 = T_0_21 * (1 + 1/eta_LPC * (pi_lpc**((kappa_air - 1) / kappa_air) - 1))

    # HPC 25 -> 3
    P_0_3 = P_0_25 * pi_hpc
    T_0_3 = T_0_25 * (1 + 1/eta_HPC * (pi_hpc**((kappa_air - 1) / kappa_air) - 1))

    # Combustor 3 -> 4
    T_04 = T_04
    m_dot_fuel = (m_dot_c * cp_gas * (T_04 - T_0_3)) / (eta_combustor * LHV)
    m_dot_c2 = m_dot_c + m_dot_fuel
    p_04 = P_0_3 * pi_comb

    # Continue
    return eta_comb, eta_th, eta_jet, eta_prop, eta_tot