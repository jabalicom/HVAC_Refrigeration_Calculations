# -*- coding: utf-8 -*-
"""

Script for  calculating Air Psychrometric Proceses for HVAC and Thermal Comfort. Air properties are calculated 
using Coolprop library. 
The script calculates the humidity ratio, air enthalpy, sensible heat, latent heat and total heat added/removed
from the air based on the input varialbles (input and output temperatures, relative humidities and airflow speed
and duct area).
The psychrometric process is them plotted on a psychrometric chart and saved as .svg file.

Additionally the script calculates:
    - Adiabatic cooling of air
    - Adiabatic mixing of air streams
    - Air inlet and outlet dew point and frost points are calculated as well.
    
Note: atmospheric pressure assumed to be 101325 Pa.

"""
import numpy as np
import matplotlib.pyplot as plt
from CoolProp.HumidAirProp import HAPropsSI
import CoolProp.CoolProp as CP
import math

# Input Variables

ti = 31 + 273.15 # Inlet temperature 
rhi = 50 /100 # Inlet Relative Humidity 
to = 22 + 273.15 # Oulet temperature 
rho = 50 /100 # Outlet Relative Humidity 
air_velocity = 2.09 # Air duct Velocity [m/s]
duct_area = 0.2 * 0.2 #[m2]

# Constants

Cp = 1.006 # Specific Heat of Air [kJ/kgK]
hwe = 2454 # Latent heat evaporization water- in air at atmospheric pressure and 20oC- [kJ/kg]
p_air = 1.202 # Air density [kg/m3]

print('Inlet Temperature = {:.2f}K or {:.2f} Celsius'.format(ti, ti-273.15))
print('Inlet %RH = {:.2f}%'.format(rhi*100))
print('Outlet Temperature = {:.2f}K or {:.2f} Celsius'.format(to, to-273.15))
print('Outlet %RH = {:.2f}%'.format(rho*100))
print('Air Velocity = {:.2f} m/s'.format(air_velocity)) 
print('Duct Area = {:.2f} m2'.format(duct_area))
print('\n')

# Heating or Cooling Calculation

def heating_and_cooling(ti, rhi, to, rho, air_velocity, duct_area):
    
    V = duct_area * air_velocity # Volumetric Flow Rate [m3/s]
    pv = HAPropsSI('P_w','T',ti,'P',101325,'R',rhi)/100000 # Partial pressure of water vapor
    ma = ((1.0135 - pv) * V)*100000 / (287 * ti) # Air Mass Flow [kg/s]

    global Wi
    global Wo
    
    Wi = HAPropsSI('W','T',ti,'P',101325,'R',rhi) # Inlet Humidity Ratio[kg/kg of dry air]
    Wo = HAPropsSI('W','T',to,'P',101325,'R',rho) # Outlet Humidity Ratio[kg/kg of dry air]

    hi = HAPropsSI('H','T',ti,'P',101325,'R',rhi)/1000 # Inlet Air Enthalpy [kJ/kg]
    ho = HAPropsSI('H','T',to,'P',101325,'R',rho)/1000 # Outlet Air Enthalpy [kJ/kg]

    Qs = Cp * ma * (to - ti) #Sensible Heat Added/Removed from the Air [kW]
    Ql = hwe * ma *(Wo - Wi) # Latent Heat Added/Removed from the Air [kW]
    Qt = ma * (ho - hi) # Total Heat Added/Removed from the Air [kW]
    
    print('Partial Pressure of Water Vapour = {:.3f} bar'.format(pv))
    print('Mass of dry air = {:.2f} kg/s'.format(ma))

    print('Inlet Air Humidity Ratio = {:.5f} kg/kg of dry air' .format(Wi))
    print('Outlet Air Humidity Ratio = {:.5f} kg/kg of dry air' .format(Wo))

    print('Inlet Air Enthalpy = {:.3f} kJ/kg of dry air' .format(hi))
    print('Outlet Air Enthalpy = {:.3f} kJ/kg of dry air' .format(ho))

    print('Sensible Heat added/removed = {:.3f} kJ/s or kW'.format(Qs))
    print('Latent Heat added/removed = {:.3f} kJ/s or kW'.format(Ql))
    print('Amount of Heat added/removed = {:.2f} kJ/s or kW'.format(Qt))
    print('\n')
    
# Psychrometric Process Plot

def plot_psychrometrics(ti, Wi, to, Wo ):
    
    fig, ax = plt.subplots(1,1,figsize=(10, 8))
    Tdbvec = np.linspace(-10, 55)+273.15

    plt.plot([ti-273.15, to-273.15], [Wi,Wo], color = 'g', lw = 5, marker = 'o', ms = 10, mfc = 'y')

    # Lines of constant relative humidity
    for RH in np.arange(0.1, 1, 0.1):
        W = CP.HAPropsSI("W","R",RH,"P",101325,"T",Tdbvec)
        plt.plot(Tdbvec-273.15, W, color='k', lw = 0.5)

    # Saturation curve
    W = CP.HAPropsSI("W","R",1,"P",101325,"T",Tdbvec)
    plt.plot(Tdbvec-273.15, W, color='k', lw=1.5)

    plt.plot([10, 20], [0.05, 0.1], color='k', lw=1.5)
    # Lines of constant Vda
    for Vda in np.arange(0.69, 0.961, 0.01):
        R = np.linspace(0,1)
        W = CP.HAPropsSI("W","R",R,"P",101325,"Vda",Vda)
        Tdb = CP.HAPropsSI("Tdb","R",R,"P",101325,"Vda",Vda)
        plt.plot(Tdb-273.15, W, color='b', lw=1.5 if abs(Vda % 0.05) < 0.001 else 0.5)

    # Lines of constant wetbulb
    for Twb_C in np.arange(-16, 33, 2):
        if Twb_C == 0:
            continue
        R = np.linspace(0.0, 1)
        #print(Twb_C)
        Tdb = CP.HAPropsSI("Tdb","R",R,"P",101325,"Twb",Twb_C+273.15)
        W = CP.HAPropsSI("W","R",R,"P",101325,"Tdb",Tdb)
        plt.plot(Tdb-273.15, W, color='r', lw=1.5 if abs(Twb_C % 10) < 0.001 else 0.5)
        plt.plot([10, 20], label="Line 1")


    plt.xlabel(r'Dry bulb temperature $T_{\rm db}$ ($^{\circ}$ C)')
    plt.ylabel(r'Humidity Ratio $W$ (kg/kg)')
    plt.ylim(0, 0.030)
    plt.xlim(-10, 55)
    plt.savefig('Psychrometric Process.svg')
    plt.show()

# Dew Point Calculation 
def get_dew_point_c(ti, rhi):

    A = 17.27
    B = 237.7
    t_air_c = ti - 273.15
    alpha = ((A * t_air_c) / (B + t_air_c)) + math.log(rhi)
    global dew_point_C
    dew_point_C = (B * alpha) / (A - alpha)
    print('Inlet Air Dew Point = {}C'.format(round(dew_point_C,2)))

# Frost Point Calculation
def get_frost_point_c(ti, dew_point_c):
 
    dew_point_k = 273.15 + dew_point_c 
    frost_point_k = dew_point_k - ti + 2671.02 / ((2954.61 / ti) + 2.193665 * math.log(ti) - 13.3448)
    global frost_point_C
    frost_point_C = frost_point_k - 273.15
    print('Inlet Air Frost Point = {:.2f}C'.format(frost_point_C))

# Adiabatic Cooling of Air:
def adiabatic_cooling(ti, rhi, to, rho):
    
    
    V = duct_area * air_velocity * 60 * 60 #Volumetric Air Flow [m3/h]
    pv = HAPropsSI('P_w','T',to,'P',101325,'R',rho)/100000 # Partial pressure of water vapor
    va = (287 * ti) /((1.0135 - pv)*100000) #[m3/kg of dry air]

    hi = HAPropsSI('H','T',ti,'P',101325,'R',rhi)/1000 # Inlet Air Enthalpy[kJ/kg of dry air]
    ho = HAPropsSI('H','T',to,'P',101325,'R',rho)/1000 # Outlet Air Enthalpy[kJ/kg of dry air]

    W1 = HAPropsSI('W','T',ti,'P',101325,'R',rhi) #[kg/kg of dry air]
    W2 = HAPropsSI('W','T',to,'P',101325,'R',rho)   #[kg/kg of dry air]

    spray_water = (W2 - W1) / va #
    amount_of_water_required_per_hour = spray_water * V #[kg/h]
    Qadiabatic = abs(p_air * (ho-hi) * V/3600)

    print('\n')
    print('Air Volumetric Rate = {:.2f} m3/h'.format(V))
    print('Inlet Air Enthalpy = {:.3f} kJ/kg of dry air' .format(hi))
    print('Outlet Air Enthalpy = {:.3f} kJ/kg of dry air' .format(ho))
    print('Partial Pressure of Water Vapour = {:.3f} bar'.format(pv))
    print('Specific Humidity Inlet Air = {:.5f} kg/kg of dry air'.format(W1))
    print('Specific Humidity Outlet Air = {:.5f} kg/kg of dry air'.format(W2))
    print('Spray Water = {:.5f} kg of moisture/m3 of air'.format(spray_water))
    print('Amount of Water Required per hour = {:.3f} kg/h'.format(amount_of_water_required_per_hour))
    print('Adiabatic Cooler Power Output = {:.2f} kW'.format(Qadiabatic))
    print('\n')
 

# Inddor Thermal Comfort PPD and PMV Calculation
    
def comfPMV(ta, tr, vel, rh, met, clo, wme):
    """
    returns [pmv, ppd]
    ta, air temperature (C)
    tr, mean radiant temperature (C)
    vel, relative air velocity (m/s)
    rh, relative humidity (%) Used only this way to input humidity level
    met, metabolic rate (met)--> Seated=1;Standing=1.3;Walking=1.7;medium Work=3
    clo, clothing (clo)--> Typical Summer = 0.5;Typical Winter = 1
    wme, external work, normally around 0 (met)
    
    PMV Scale--> -3=Cold;-2=Cool;-1=Slightly Cool;0=Neutral;1=Slightly Warm;2=Warm;3=Hot
    5% < PPD < 100%
    
    """

    pa = rh * 10 * math.exp(16.6536 - 4030.183 / (ta + 235))

    icl = 0.155 * clo  # thermal insulation of the clothing in M2K/W
    m = met * 58.15  # metabolic rate in W/M2
    w = wme * 58.15  # external work in W/M2
    mw = m - w  # internal heat production in the human body
    if (icl <= 0.078):
        fcl = 1 + (1.29 * icl)
    else:
        fcl = 1.05 + (0.645 * icl)

    # heat transf. coeff. by forced convection
    hcf = 12.1 * math.sqrt(vel)
    taa = ta + 273
    tra = tr + 273
    tcla = taa + (35.5 - ta) / (3.5 * icl + 0.1)

    p1 = icl * fcl
    p2 = p1 * 3.96
    p3 = p1 * 100
    p4 = p1 * taa
    p5 = (308.7 - 0.028 * mw) + (p2 * math.pow(tra / 100, 4))
    xn = tcla / 100
    xf = tcla / 50
    eps = 0.00015

    n = 0
    while abs(xn - xf) > eps:
        xf = (xf + xn) / 2
        hcn = 2.38 * math.pow(abs(100.0 * xf - taa), 0.25)
        if (hcf > hcn):
            hc = hcf
        else:
            hc = hcn
        xn = (p5 + p4 * hc - p2 * math.pow(xf, 4)) / (100 + p3 * hc)
        n += 1
        if (n > 150):
            print ('Max iterations exceeded')
            return 1


    tcl = 100 * xn - 273

    # heat loss diff. through skin
    hl1 = 3.05 * 0.001 * (5733 - (6.99 * mw) - pa)
    # heat loss by sweating
    if mw > 58.15:
        hl2 = 0.42 * (mw - 58.15)
    else:
        hl2 = 0
    # latent respiration heat loss
    hl3 = 1.7 * 0.00001 * m * (5867 - pa)
    # dry respiration heat loss
    hl4 = 0.0014 * m * (34 - ta)
    # heat loss by radiation
    hl5 = 3.96 * fcl * (math.pow(xn, 4) - math.pow(tra / 100, 4))
    # heat loss by convection
    hl6 = fcl * hc * (tcl - ta)

    ts = 0.303 * math.exp(-0.036 * m) + 0.028
    pmv = ts * (mw - hl1 - hl2 - hl3 - hl4 - hl5 - hl6)
    ppd = 100.0 - 95.0 * math.exp(-0.03353 * pow(pmv, 4.0)
        - 0.2179 * pow(pmv, 2.0))

    r = []
    r.append(pmv)
    r.append(ppd)

    print('Thermal Comfort Predicted Mean Vote = {:.2f}'.format(pmv))
    print('Thermal Comfort Percentage of People Disatisfied = {:.2f}%'.format(ppd))

    
if __name__ == "__main__": 
    
    heating_and_cooling(ti, rhi, to, rho, air_velocity, duct_area)
    get_dew_point_c(ti, rhi)
    get_frost_point_c(ti, dew_point_C)
    plot_psychrometrics(ti, Wi, to, Wo)
    adiabatic_cooling(ti, rhi, to, rho)
    comfPMV(ta=ti-273.15, tr=ti-273.15, vel=0.1, rh=50, met=1.2, clo=0.75, wme=0)

