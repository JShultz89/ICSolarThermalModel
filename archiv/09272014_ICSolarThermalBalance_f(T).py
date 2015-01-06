# Written by Justin Shultz from Rensselaer Polytechnic Institute 
# for Assad Oberai, with his assistance too

import numpy as np

m = 3 # Number of Modules
A = 1.0 # Area m^2

q_receiver = 200.0 # Heat flow into water from Module Heat Receiver
q_module = 50.0 # Heat flow into air from Heat Loss from the Module

C_pw = 4.186 # Specific heat water
C_pa = 1.01 # Specific heat capacity of air

m_w = 2.0 # Mass flowrate of water
m_a = 4.0 # Flowrate of air

h_a = 3.66 # Convection heat transfer coefficient between water and interior of pipe
h_wa = 3.66 # Convection heat transfer coefficient between air and water pipe
h_i = 3.66 # Convection heat transfer coefficient between air and interior
h_e = 3.66 # Convection heat transfer coefficient between air and exterior

T_i = 20.0 # Temperature on the interior of the building, degrees C
T_e = 25.0 # Temperature on the exterior of the building, degrees C
T_ai = 20.0 # Inlet temperature of air, degrees C
T_wi = 25.0 # Inlet temperature of water, degrees C
T = [] # Define T as an array
balance = []

# This function creates an array of temperature values (Even is air, Odd is water)
# The values are an initial guess for use in the residue function
def GenTemperatureArray(m):
    count = 0
    T = []
    T_a_tmp = T_ai # Temperary value for the water
    T_w_tmp = T_wi # Temperary value for the water 
    guessIncrement = 5 # Increases temperature in each stage by number, may not be necessary
    while count < 4*m+1:
        if count%2 == 0:
            T_tmp = T_a_tmp
            T.append(T_tmp) # Add the first air temperature to the array
            T_a_tmp += guessIncrement # Increase each air temperature by 5 degrees 
        else:
            T_tmp = T_w_tmp
            T.append(T_tmp) # Add the first water temeperature to the end of the array
            T_w_tmp += guessIncrement # Increase each water temperature by 5 degrees
        count += 1
    return T
    #    print T

# Finding solution in which i (for water) and i+1 (for air) returns a 0 result
#def residue
def Residue(m,T):
    for i in range(0,m):
        print "Running module %d" %(i+1)
        # Odd region (transfer between water, air, interior, and exterior)
        # Calculate water temperature for region 1 (Transfer with air)
        water_r1 = m_w * C_pw * (T_w[i+1]-T_w[i]) - h_wa*A*((T_a[i+1]-T_a[i]/2.0-(T_w[i+1]-T_w[i])/2.0))
        # Calculate air temperature for region 1 (Transfer with water, interior and exterior)
        air_r1 = m_a * C_pa*(T_a[i+1]-T_a[i]) + h_wa*A*((T_a[i+1]-T_a[i])/2.0-(T_w[i+1]-T_w[i])/2.0)-h_i*A*(T_i-(T_a[i+1]-T_a[i])/2.0)-h_e*A*(T_e-(T_a[i+1]-T_a[i])/2.0)
        
        # Even Region (Heat input)
        water_r2 = m_w * C_pw * (T_w[i+1]-T_w[i]) + q_receiver
        air_r2 = m_a * C_pa*(T_a[i+1]-T_a[i]) + q_module
        
        result = water_r1 + air_r1 + water_r2 + air_r2 # This equation should equal zero
        balance.append(result)
        
print GenTemperatureArray(3)