# This file calculates the efficiencies of a turbofan engine to then plot a sankey diagram
# Authors: Alfonso Rato, Jakub Bartmański
import inspect
import numpy as np
import matplotlib.pyplot as plt
from typing import Literal
from matplotlib.sankey import Sankey

class GasTurbine(object):
   
    # TODO: CHECK IF THE UNITS ARE CONSISTENT

    # CONSTANTS

    # Reference freestream static conditions (ISA)
    T_0_ISA = 288.15    # K
    P_0_ISA = 101325    # Pa
    R_AIR = 287         # J/(kg*K)

    def __init__(self, M_0, T_0, P_0, *, cp_air=1000, kappa_air=1.4, cp_gas=1150, kappa_gas=1.33):

        """
        Initialize the GasTurbine with given parameters.

        Parameters:
        -----------
        - M_0 (float): Freestream Mach number.
        - T_0 (float): Freesteam static temperature of the ambient air.
        - P_0 (float): Freestream static pressure of the ambient air.

        Keyword Arguments:
        ------------------
        - cp_air (float, optional): Specific heat capacity of air at constant pressure. Default is 1000 J/(kg*K).
        - kappa_air (float, optional): Adiabatic index of air. Default is 1.4.
        - cp_gas (float, optional): Specific heat capacity of gas at constant pressure. Default is 1150 J/(kg*K).
        - kappa_gas (float, optional): Adiabatic index of gas. Default is 1.33.
        """
        
        self.M_0 = M_0
        self.T_0 = T_0
        self.P_0 = P_0
        self.cp_air = cp_air
        self.kappa_air = kappa_air
        self.cp_gas = cp_gas
        self.kappa_gas = kappa_gas

        # Resolved flag
        self.__is_resolved = False

        # Efficiencies
        self.__efficiencies = None
        
    def resolve(self, \
                bpr, \
                m_dot_core_isa, \
                pi_fan, pi_hpc, pi_comb, opr, \
                eta_mech, eta_isen_intake, eta_isen_comp, eta_combustor, eta_isen_turb, eta_isen_nozzle, \
                T_0_4, LHV, *, \
                verbose=True):
        
        # Assert that the GasTurbine object is not resolved
        assert not self.__is_resolved, "GasTurbine object is already resolved. Create a new object or reset the current one."

        # Retrieve the constants and attributes
        T_0_ISA = GasTurbine.T_0_ISA
        P_0_ISA = GasTurbine.P_0_ISA
        R_AIR = GasTurbine.R_AIR

        cp_air = self.cp_air 
        kappa_air = self.kappa_air
        cp_gas = self.cp_gas
        kappa_gas = self.kappa_gas

        M_0 = self.M_0
        T_0 = self.T_0
        P_0 = self.P_0

        # Begin the calculations

        # Calculate pressure ratio acroos LPC from Overall Pressure Ratio (OPR)
        pi_lpc = opr / (pi_fan * pi_hpc)

        # Calculate freestream velocity
        V_0 = M_0 * (kappa_air * R_AIR * T_0) ** 0.5

        # Calculate core massflow from corrected core massflow
        m_dot_core = m_dot_core_isa * P_0 / P_0_ISA * (T_0_ISA / T_0) ** 0.5

        # Calculate bypass massflow from bypass ratio
        m_dot_bp = bpr * m_dot_core

        # Calculate total temperature and pressure at freestream (0 station is the freestream)
        T_0_0 = T_0 * (1 + (kappa_air - 1) / 2 * M_0 ** 2)
        P_0_0 = P_0 * (1 + (kappa_air - 1) / 2 * M_0 ** 2) ** (kappa_air / (kappa_air - 1))

        # print(f"--------------\nFreestream conditions:\nT_00: {T_0_0:.2f} K\nP_00: {P_0_0:.2f} Pa\n")

        # Inlet 0 -> 2
        P_0_2 = P_0 * (1 + eta_isen_intake * (kappa_air - 1) / 2 * M_0 ** 2) ** (kappa_air / (kappa_air - 1))
        T_0_2 = T_0_0

        # print(f"--------------\nInlet conditions:\nT_02: {T_0_2:.2f} K\nP_02: {P_0_2:.2f} Pa\n")

        # Fan 2 -> 21
        P_0_21 = P_0_2 * pi_fan
        T_0_21 = T_0_2 * (1 + 1 / eta_isen_comp * (pi_fan ** ((kappa_air - 1) / kappa_air) - 1))

        # print(f"--------------\nFan conditions:\nT_021: {T_0_21:.2f} K\nP_021: {P_0_21:.2f} Pa\n")

        # Work of fan
        W_dot_fan = (m_dot_core + m_dot_bp) * cp_air * (T_0_21 - T_0_2)

        # Bypass 2 -> 13
        P_0_13 = P_0_21
        T_0_13 = T_0_21

        # LPC 21 -> 25
        P_0_25 = P_0_21 * pi_lpc
        T_0_25 = T_0_21 * (1 + 1 / eta_isen_comp * (pi_lpc ** ((kappa_air - 1) / kappa_air) - 1))

        # print(f"--------------\nLPC conditions:\nT_025: {T_0_25:.2f} K\nP_025: {P_0_25:.2f} Pa\n")

        # Work of LPC
        W_dot_LPC = m_dot_core * cp_air * (T_0_25 - T_0_21)

        # HPC 25 -> 3
        P_0_3 = P_0_25 * pi_hpc
        T_0_3 = T_0_25 * (1 + 1 / eta_isen_comp * (pi_hpc ** ((kappa_air - 1) / kappa_air) - 1))

        # print(f"--------------\nHPC conditions:\nT_03: {T_0_3:.2f} K\nP_03: {P_0_3:.2f} Pa\n")

        # Work of HPC
        W_dot_HPC = m_dot_core * cp_air * (T_0_3 - T_0_25)

        # Combustor 3 -> 4
        # The stagnation temperature at the outlet of CC: T_0_4 which is given
        m_dot_fuel = (m_dot_core * cp_gas * (T_0_4 - T_0_3)) / (eta_combustor * LHV)
        m_dot_core_tot = m_dot_core + m_dot_fuel
        P_0_4 = P_0_3 * pi_comb

        # print(f"--------------\nCombustor conditions:\nT_04: {T_0_4:.2f} K\nP_04: {P_0_4:.2f} Pa\n")

        # HPT 4 -> 44

        # Work of HPT
        W_dot_HPT = W_dot_HPC / eta_mech

        # Rest of HPT
        T_0_44 = T_0_4 - W_dot_HPT / (m_dot_core_tot * cp_gas)
        P_0_44 = P_0_4 * (1 - (1 / eta_isen_turb) * (1 - T_0_44 / T_0_4)) ** (kappa_gas / (kappa_gas - 1))

        # print(f"--------------\nHPT conditions:\nT_044: {T_0_44:.2f} K\nP_044: {P_0_44:.2f} Pa\n")

        # LPT 44 -> 5

        # Work of LPT
        W_dot_LPT = (W_dot_LPC + W_dot_fan) / eta_mech

        # Rest of LPT
        T_0_5 = T_0_44 - W_dot_LPT / (m_dot_core_tot * cp_gas)
        P_0_5 = P_0_44 * (1 - (1 / eta_isen_turb) * (1 - T_0_5 / T_0_44)) ** (kappa_gas / (kappa_gas - 1))

        # The core nozzle. First, check if the nozzle is choked.

        # Critical pressure ratio
        pi_nozzle_crit = (1 - (1 / eta_isen_nozzle) * ((kappa_gas - 1) / (kappa_gas + 1))) ** (-kappa_gas / (kappa_gas - 1))

        # Real pressure ratio
        pi_nozzle_real = P_0_5 / P_0

        # Assign a flag if the nozzle is choked
        self.is_nozzle_choked = pi_nozzle_real > pi_nozzle_crit
        is_nozzle_choked = self.is_nozzle_choked
        
        # Exit of the nozzle

        # Choked case
        if is_nozzle_choked:

            if verbose:
                print("Core nozzle is choked.")

            # Note these are static values! 
            P_8 = P_0_5 / pi_nozzle_crit
            T_8 = T_0_5 * 2 / (kappa_gas + 1)
            rho_8 = P_8 / (R_AIR * T_8)
            V_8 = (kappa_gas * R_AIR * T_8) ** 0.5
            A_8 = m_dot_core_tot / (rho_8 * V_8)

        # Unchoked case
        else:

            if verbose:
                print("Core nozzle is unchoked.")
        
            P_8 = P_0
            T_8 = T_0_5 * (1 - eta_isen_nozzle * (1 - (P_0 / P_0_5) ** ((kappa_gas - 1) / kappa_gas)))
            rho_8 = P_8 / (R_AIR * T_8)
            V_8 = (2 * cp_gas * (T_0_5 - T_8)) ** 0.5
            A_8 = m_dot_core_tot / (rho_8 * V_8)

        # Jet effective velocity of the nozzle
        V_8_eff = V_8 + (P_8 - P_0) / (rho_8 * V_8)

        # Now do the same for bypass (nozzle) (station 18)

        # Critical pressure ratio
        pi_bp_crit = (1 - (1 / eta_isen_nozzle) * ((kappa_air - 1) / (kappa_air + 1))) ** (-kappa_air / (kappa_air - 1))

        # Real pressure ratio
        pi_bp_real = P_0_13 / P_0

        # Assign a flag if the bypass nozzle is choked
        self.is_bp_choked = pi_bp_real > pi_bp_crit
        is_bp_choked = self.is_bp_choked
        
        # Exit of the bypass nozzle

        # Choked case
        if is_bp_choked:

            if verbose:
                print("Bypass nozzle is choked.")
    
            P_18 = P_0_13 / pi_bp_crit
            T_18 = T_0_13 * 2 / (kappa_air + 1)
            rho_18 = P_18 / (R_AIR * T_18)
            V_18 = (kappa_air * R_AIR * T_18) ** 0.5
            A_18 = m_dot_bp / (rho_18 * V_18)

        # Unchoked case
        else:

            if verbose:
                print("Bypass nozzle is unchoked.")

            P_18 = P_0
            T_18 = T_0_13 * (1 - eta_isen_nozzle * (1 - (P_0 / P_0_13) ** ((kappa_air - 1) / kappa_air)))
            rho_18 = P_18 / (R_AIR * T_18)
            V_18 = (2 * cp_air * (T_0_13 - T_18)) ** 0.5
            A_18 = m_dot_bp / (rho_18 * V_18)

        # Jet effective velocity of the bypass nozzle
        V_18_eff = V_18 + (P_18 - P_0) / (rho_18 * V_18)

        # Calculate gas generator power
        W_dot_fan_core = m_dot_core * cp_air * (T_0_21 - T_0_2)
        T_0_g = (W_dot_fan_core + W_dot_LPC) / (m_dot_core_tot * cp_gas * eta_mech)
        P_0_g = P_0_44 * (1 - ((1 - T_0_g / T_0_44) / eta_isen_turb)) ** (kappa_gas / (kappa_gas - 1))
        W_dot_gg = m_dot_core_tot * cp_gas * T_0_g * (1 - (P_0 / P_0_g) ** ((1- kappa_gas) / kappa_gas))# TODO

        # Finally, calculate the efficiencies

        # Combustion efficiency is given
        eta_comb = eta_combustor

        # Thermodynamic efficiency
        eta_td = W_dot_gg / (m_dot_core * cp_gas * (T_0_4 - T_0_3))

        # Jet efficiency
        eta_jet = (0.5 * m_dot_bp * (V_18_eff ** 2 - V_0 ** 2) + 0.5 * m_dot_core_tot * (V_8_eff ** 2 - V_0 ** 2)) / W_dot_gg

        # Propulsive efficiency
        eta_prop = V_0 * (m_dot_bp * (V_18_eff - V_0) + m_dot_core_tot * (V_8_eff - V_0)) \
                    / (0.5 * m_dot_bp * (V_18_eff ** 2 - V_0 ** 2) + 0.5 * m_dot_core_tot * (V_8_eff ** 2 - V_0 ** 2))
        
        # Thermal efficiency
        eta_th = (0.5 * m_dot_bp * (V_18_eff ** 2 - V_0 ** 2) + 0.5 * m_dot_core_tot * (V_8_eff ** 2 - V_0 ** 2)) / (m_dot_fuel * LHV)

        # Total efficiency
        eta_tot = eta_th * eta_prop

        # Store the efficiencies
        self.__efficiencies = eta_comb, eta_td, eta_jet, eta_prop, eta_th, eta_tot

        print(f"--------------\nCombustion efficiency: {eta_comb:.2f}\n")
        print(f"Thermodynamic efficiency: {eta_td:.2f}")
        print(f"Jet efficiency: {eta_jet:.2f}\n")
        print(f"Propulsive efficiency: {eta_prop:.2f}")
        print(f"Thermal efficiency: {eta_th:.2f}\n")
        print(f"Total efficiency: {eta_tot:.2f}")


        # Set the resolved flag
        self.__is_resolved = True
        
        return eta_comb, eta_td, eta_jet, eta_prop, eta_th, eta_tot
    
    def reset(self):     
            # Reset the resolved flag
            self.__is_resolved = False
    
            # Reset the efficiencies
            self.__efficiencies = None

    def draw_sankey(self):
        
        # Assert that the GasTurbine object is resolved
        assert self.__is_resolved, "GasTurbine object must be resolved before drawing the sankey diagram. Call the resolve method first."

        # Retrieve efficiencies values
        eta_comb, eta_td, eta_jet, eta_prop, eta_th, eta_tot = self.__efficiencies

        # Pack the values into numpy array in the order of the sankey diagram:
        # eta_comb, eta_td, eta_jet, eta_prop
        values = np.array([eta_comb, eta_td, eta_jet, eta_prop]) * 100

        # Invert the values for the losses
        values = 100 - values

        # print(values)

        fig = plt.figure(figsize=(8.3, 6))
        ax = fig.add_subplot(1, 1, 1, xticks=[], yticks=[])
        
        # Define the flows and layout
        sankey = Sankey(ax=ax, scale=0.011, head_angle=90, format='%.0f', unit='%')

        # Add the flows with improved offsets and black outlines
        sankey.add(flows=[100, -values[0], -values[1], -values[2], -values[3], values.sum() - 100],  # Positive = forward flow, negative = losses
                labels=[None, None, None, None, None, None],  # Labels handled manually
                orientations=[0, 1, 1, 1, 1, 0],  # -1 = downward, 0 = rightward
                pathlengths=[0.5, 1., 1., 1., 1., 1.],  # Adjust lengths for separation
                facecolor="#ffcccc",  # Light pink
                edgecolor="black")   # Pure black outline

        # Render the Sankey diagram
        diagrams = sankey.finish()

        # Construct labels from values
        l1 = f"{values[0]:.1f}%"
        l2 = f"{values[1]:.1f}%"
        l3 = f"{values[2]:.1f}%"
        l4 = f"{values[3]:.1f}%"
        l5 = f"{100-values.sum():.1f}%"

        # Manually add text in boxes
        labels = [
            ("Chemical Power", "100%"),
            ("Incomplete\n Combustion", l1),
            ("Heat Loss\n (Heat)", l2),
            ("Heat Loss\n (Gas Power)", l3),
            ("KE Loss", l4),
            ("Thrust Power", l5)
        ]

        y_pos = 0.7
        positions = [
            (-0.69, 0.),     # Chemical Power
            (0.4, y_pos),     # Incomplete Combustion
            (1., y_pos),   # Heat Loss (Heat)
            (1.4, y_pos),   # Heat Loss (Gas Power)
            (1.8, y_pos),   # Kinetic Energy Loss
            (2.6, -0.33)    # Thrust Power
        ]

        for (label, value), (x, y) in zip(labels, positions):
            ax.text(
                x, y, f"{label}\n{value}", 
                bbox=dict(facecolor="white", edgecolor="black", boxstyle="round,pad=0.3"),
                ha="center", va="center", fontsize=9
            )

        ax.title.set_fontsize(14)

        plt.tight_layout()
        plt.show()

class GasTurbineData(object):

    def __init__(self):
        
        # Common characteristics
        self.common = {
            'LHV': 43E6,    # J/kg
            'eta_isen_intake': 0.99,    # -
            'eta_mech': 0.99,    # -
            'pi_comb': 0.96,    # -
            'eta_isen_nozzle': 0.99,    # -
        }

        # JT8D characteristics
        self.jt8d = {
            'bpr': 1.62,    # -
            'm_dot_core_isa': 90.2,    # kg/s
            'pi_fan': 1.9,    # -
            'pi_hpc': 3.5,    # -
            'T_0_4': 1150,    # K
            'opr': 17,    # -
            'eta_isen_comp': 0.85,    # -
            'eta_isen_turb': 0.88,    # -
            'eta_combustor': 0.985    # -
        } # TODO: WE'RE NOT USING STAGE INFO OR ENGINE DIAMETER?

        # LEAP-1B characteristics
        self.leap1b = {
            'bpr': 8.6,    # -
            'm_dot_core_isa': 50.0,    # kg/s
            'pi_fan': 1.5,    # -
            'pi_hpc': 10.0,    # -
            'T_0_4': 1450,    # K
            'opr': 40,    # -
            'eta_isen_comp': 0.92,    # -
            'eta_isen_turb': 0.92,    # -
            'eta_combustor': 0.995    # -
        } # TODO: WE'RE NOT USING STAGE INFO OR ENGINE DIAMETER?

        # # Check if jt8d dictionary keys are all the same as position arguments of GasTurbine.resolve
        # for key in inspect.signature(GasTurbine.resolve).parameters:
        #     print(f"Key '{key}' not found in jt8d dictionary") if key not in jt8d.keys() else None

    def get(self, type: Literal['jt8d', 'leap1b']):
        
        # Assert that the type is either 'jt8d' or 'leap1b'
        assert type in ['jt8d', 'leap1b'], f"Type '{type}' not found in GasTurbineData"

        # Get the main characteristics
        if type == 'jt8d':
            main = self.jt8d
        elif type == 'leap1b':
            main = self.leap1b
        # ...

        # Returned merged dictionary
        return main | self.common

if __name__ == """__main__""":

    # Ambient conditions
    ambient = {
        'T_0': 220,     # K
        'M_0': 0.78,    # Mach number
        'P_0': 23842    # Pa
    } # TODO: WE'RE NOT USING ALTIDUTE?

    # type = 'leap1b'
    type = 'jt8d'

    # Get data for appropriate engine type
    data = GasTurbineData().get(type)

    # Create the GasTurbine object and resolve it
    gt = GasTurbine(**ambient)
    gt.resolve(**data)
    # gt.draw_sankey()