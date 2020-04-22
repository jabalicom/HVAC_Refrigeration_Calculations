# -*- coding: utf-8 -*-
"""
Function to evaluate performance of a refrigeration cycle based on readings
that can/commonly taken "on site"
The script also plots the P-h and T-s cycle of the refrigeration unit.
The script uses CoolProp lybrary for thermodynamic calculations. 
Note: atmospheric pressure assumed to be 101325 Pa
"""

import CoolProp.CoolProp as CP
from CoolProp.Plots import PropertyPlot
from CoolProp.Plots import SimpleCompressionCycle
import matplotlib.pyplot as plt

# Refrigeration Input Variables

fluid = 'R134a'
T_suction = -0.15 
T_discharge = 69.2 
T_TEVinlet = 31.5 
T_TEVoutlet = -2.6
P_low_mano = 1.3
P_high_mano = 7.94

def refrigeration_cycle(fluid,T_suction,T_discharge,T_TEVinlet,T_TEVoutlet,P_low_mano,P_high_mano):
    
    T_K_suction = T_suction +273.15
    T_K_discharge = T_discharge + 273.15
    T_K_TEVinlet = T_TEVinlet + 273.15
    #T_K_TEVoutlet = T_TEVoutlet + 273.15

    P_sat_evaporating = (P_low_mano + 1) * 100000
    P_sat_condensing = (P_high_mano + 1) * 100000
    P_ratio = P_sat_condensing / P_sat_evaporating 

    T_K_sat_evaporating = CP.PropsSI('T','P',P_sat_evaporating,'Q',0,fluid)
    T_K_sat_condensing = CP.PropsSI('T','P',P_sat_condensing,'Q',1,fluid)

    superheat = T_K_suction - T_K_sat_evaporating
    subcooling = T_K_sat_condensing - T_K_TEVinlet

    h_suction = CP.PropsSI('H','P',P_sat_evaporating,'T',T_K_suction,fluid)/1000
    h_discharge = CP.PropsSI('H','P',P_sat_condensing,'T',T_K_discharge,fluid)/1000
    h_TEVinlet = CP.PropsSI('H','P',P_sat_condensing,'T',T_K_TEVinlet,fluid)/1000


    density_liquid_refrigerant = CP.PropsSI('D','P',P_sat_condensing,'T',T_K_TEVinlet,fluid)/1000
    Q_refrigeration = h_suction - h_TEVinlet
    Q_compression = h_discharge - h_suction
    Q_condenser = h_TEVinlet - h_discharge

    s_suction = CP.PropsSI('S','P',P_sat_evaporating,'T',T_K_suction,fluid)/1000 
    h_isentropic = CP.PropsSI('H','P',P_sat_condensing,'S',s_suction*1000,fluid)/1000 
    h_Isen_Comp_Work_Enthalpy = h_isentropic - h_suction
    #global Isentropic_Efficiency
    Isentropic_Efficiency = h_Isen_Comp_Work_Enthalpy*100 / Q_compression


# Plot P-h Cycle
    def plot_Ph():
        pp = PropertyPlot('HEOS::R134a', 'PH', unit_system='EUR', tp_limits='ACHP')
        pp.calc_isolines()
        cycle = SimpleCompressionCycle('HEOS::R134a', 'PH', unit_system='EUR', tp_limits='ACHP')
        cycle.simple_solve_dt(Te=T_K_sat_evaporating ,Tc=T_K_sat_condensing, dT_sh=superheat, dT_sc=subcooling, eta_com= Isentropic_Efficiency/100, SI=True)
        cycle.steps = 50
        sc = cycle.get_state_changes()
        plt.close(cycle.figure)
        pp.draw_process(sc)
        pp.savefig('P-h.svg')
        pp.show()


    # Plot T-s Cycle
    def plot_Ts():
        pp = PropertyPlot('HEOS::R134a', 'TS', unit_system='EUR', tp_limits='ACHP')
        pp.calc_isolines()
        cycle = SimpleCompressionCycle('HEOS::R134a', 'TS', unit_system='EUR', tp_limits='ACHP')
        cycle.simple_solve_dt(Te=T_K_sat_evaporating ,Tc=T_K_sat_condensing, dT_sh=superheat, dT_sc=subcooling, eta_com= Isentropic_Efficiency/100, SI=True)
        cycle.steps = 50
        sc = cycle.get_state_changes()
        plt.close(cycle.figure)
        pp.draw_process(sc)
        pp.savefig('T-s.svg')
        pp.show()

    plot_Ph()
    plot_Ts()
    print("Refrigerant Used = {}".format(fluid))
    print("Suction Line Temperature = {} C".format(T_suction))
    print("Discharge Line Temperature = {} C".format(T_discharge))
    print("TEV Inlet Temperature = {} C".format(T_TEVinlet))
    print("TEV Outlet Temperature = {} C".format(T_TEVoutlet))
    print("Manometric Pressure Low Side = {} bar".format(P_low_mano))
    print("Manometric Pressure High Side = {} bar".format(P_high_mano))
    print('\n')
    print("Saturated Evaporating Pressure = {} Pa".format(round(P_sat_evaporating,0)))
    print("Saturated Condensing Pressure = {} Pa".format(round(P_sat_condensing,0)))
    print("Pressure Ratio = {}".format(round(P_ratio,2)))
    print("Saturated Evaporating Temperature = {} C".format(round(T_K_sat_evaporating -273.15,2)))
    print("Saturated Condensing Temperature = {} C".format(round(T_K_sat_condensing -273.15,2)))
    print("Superheat = {} C".format(round(superheat,2)))     
    print("Subcooling = {} C".format(round(subcooling,2)))      
    print('\n')
    print("Suction Line Enthalpy = {} kJ/kg".format(round(h_suction,2)))      
    print("Discharge Line Enthalpy = {} kJ/kg".format(round(h_discharge,2)))      
    print("TEV Inlet Enthalpy = {} kJ/kg".format(round(h_TEVinlet,2)))      
    print('\n')
    print("Refrigerant Density (liquid) = {} kg/l".format(round(density_liquid_refrigerant,3)))      
    print("Refrigeration Effect per kg of refrigerant = {} kJ/kg".format(round(Q_refrigeration ,3)))      
    print("Heat of Compression per kg of refrigerant = {} kJ/kg".format(round(Q_compression,3)))
    print("Heat Rejected at Condenser per kg of refrigerant = {} kJ/kg".format(round(Q_condenser,3)))
    print('\n')
    print('Suction Line Entropy = {} kJ/kgK'.format(round(s_suction,3)))      
    print("Isentropic Enthalpy = {} kJ/kg".format(round(h_isentropic,3)))      
    print("Isentropic Compression Work Enthalpy = {} kJ/kg".format(round(h_Isen_Comp_Work_Enthalpy ,3)))      
    print("Isentropic Efficiency = {}%".format(round(Isentropic_Efficiency,3))) 


if __name__ == "__main__": 
    refrigeration_cycle(fluid,T_suction,T_discharge,T_TEVinlet,T_TEVoutlet,P_low_mano,P_high_mano)

 