# HVAC_Refrigeration_Calculations
Two Python scrips to evaluate the performance of a refrigeration cycle and calculating Air Psychrometric Proceses.
## Script: refrigeration.py
- Function to evaluate performance of a refrigeration cycle based on readings(script inputs) that are/can be commonly taken "on site"
- The script also plots the P-h and T-s cycle of the refrigeration unit.
- The script uses [CoolProp](https://github.com/CoolProp/CoolProp) for thermodynamic calculations. 
- Note: atmospheric pressure assumed to be 101325 Pa. If other, calls to CoolProp need to be change accordingly.
## Script: psychrometrics.py
- Script for  calculating Air Psychrometric Proceses for HVAC and Thermal Comfort. 
- Air properties are calculated using [CoolProp](https://github.com/CoolProp/CoolProp). 
- The script calculates the humidity ratio, air enthalpy, sensible heat, latent heat and total heat added/removed from the air based on the input varialbles (input and output temperatures, relative humidities and airflow speed and duct area).
- The psychrometric process is them plotted on a psychrometric chart and saved as .svg file.

- Additionally, the script calculates:
    - Adiabatic cooling of air
    - Adiabatic mixing of air streams
    - Air inlet and outlet dew point and frost points are calculated as well.
    
- Note: atmospheric pressure assumed to be 101325 Pa. If other, calls to CoolProp need to be change accordingly.
