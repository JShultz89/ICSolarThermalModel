# Written by Justin Shultz from Rensselaer Polytechnic Institute 
# for Assad Oberai, with his assistance too

import numpy 

m = 2 # Number of Modules
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
T_wi = 13.0 # Inlet temperature of water, degrees C

T = numpy.empty(m*4+2, dtype=float) # Define T as an array of a specifc length


def Cp(liquid, T):
    # This function calculates the Specific Heat of Water and Air based on Temperature
    # The equation for Water was a 6 degree polynomial calculated using a temperature range from 5-105C
    # The specifc heat of air does not change over the temperature range of the system and can be set constant
    if liquid == "w":
        # Constants determined from a 6 degree fit polynomial using MatLab CurveFit Toolbox
        # R = 0.9989 & SSE = 4.494e-06
        p1 = 0.0005544
        p2 = -0.001201 
        p3 = 6.553e-05
        p4 = 0.001323
        p5 = 0.008304
        p6 = 0.01152
        p7 = 4.183
        Cp = p1*T**6 + p2*T**5 + p3*T**4 + p4*T**3 + p5*T**2 + p6*T + p7
    elif liquid == "a":
        # The specifc heat of air does not change over the temperature region of the system
        Cp = 1.005
    return Cp # Return the specific heat



def h(interface,T):
    # This function returns the convection heat transfer coefficient for specific temperatures
    # Interface, at which exchange is occuring. Options: p (water interior of pipe), pa (pipe to air), i (interior), e (exterior)
    # L is the characteristic length, this changes based on which part of the system is being calculated
    
    def k(liquid, T):
        if liquid == "w":
            # Equations and values acquired from "Standard Reference Data for the Thermal Conductivity of Water" by M.L.V. Ramires 1995
            # Below equations accurate in temperature ranges from 274K to 370K
            kAtRoom = 0.6065
            dimensionless_k = -1.48445 + 4.12292*((T+274.15)/298.15) - 1.63866*((T+274.15)/298.15)**2
            k = dimensionless_k * kAtRoom
        elif liquid == "a":
            # Constants determined from a 3 degree fit polynomial using MatLab CurveFit Toolbox
            # R = 1
            p1 =  1.668e-22
            p2 =  -1.624e-20
            p3 =  7e-05
            p4 =  0.0243
            k = p1*T**3 + p2*T**2 + p3*T + p4
        return k # Return the specific heat
    
    if interface == 'p':
        k = k('w', T)
        L = 0.3
        Nu = 3.66 # for Reynold numbers < 2300 (laminar) and constant wall temperature
    elif interface == 'pa':
        k = k('a', T)
        L = 0.3
        Nu = 3.66
    elif interface == 'i':
        k = k('a', T)
        L = 3.0
        Nu = 3.66
    elif interface == 'e':
        k = k('a', T)
        L = 3.0
        Nu = 3.66
    h = k / L * Nu
    return h


def GenTemperatureArray(m):
    # This function creates an array of temperature values (Even is air, Odd is water)
    # The values are an initial guess for use in the residue function
    count = 0
    T = numpy.empty(m*4+2, dtype=float)
    T_a_tmp = T_ai # Temperary value for the water
    T_w_tmp = T_wi # Temperary value for the water 
    guessIncrement = 5 # Increases temperature in each stage by number, may not be necessary
    while count < 4*m+2:
        if count%2 == 0:
            T_tmp = T_w_tmp
            T.append(T_tmp) # Add the first air temperature to the array
            T_w_tmp += guessIncrement # Increase each air temperature by 5 degrees 
        else:
            T_tmp = T_a_tmp
            T.append(T_tmp) # Add the first water temeperature to the end of the array
            T_a_tmp += guessIncrement # Increase each water temperature by 5 degrees
        count += 1
    return T
    #    print T


def Residue(m,T):
    # Calculates the temperature balance of each region of the system
    # Inputs are temperature and number of modules
    # Calculated values are the inbalance in the heat equations
    # Inbalance will be solved using Newton's Iteration in a further step
    q = numpy.empty(m*4, dtype=float)
    j = 0
    for i in range(0,m):
        print "Running module %d" %(i+1)
        # Odd region (transfer between water, air, interior, and exterior)
        # Calculate water temperature for region 1 (Transfer with air)
        q_w1 = m_w * Cp('w',(T[j+2]-T[j])/2) * (T[j+2]-T[j]) - h('p',(T[j+3]-T[j+1])/2.0-(T[j+2]-T[j])/2) * A * ((T[j+3]-T[j+1])/2.0-(T[j+2]-T[j])/2.0)
        q[j] = q_w1
        # Calculate air temperature for region 1 (Transfer with water, interior and exterior)
        q_a1 = m_a * Cp('a',(T[j+3]-T[j+1])/2) * (T[j+3]-T[j+1]) + h('pa',(T[j+3]-T[j+1])/2.0-(T[j+3]-T[j+1])/2.0) * A * ((T[j+3]-T[j+1])/2.0 - (T[j+2]-T[j])/2.0) - h('i',T_i-(T[j+3]-T[j+1])/2) * A * (T_i-(T[j+3]-T[j+1])/2.0) - h('e',T_e-(T[j+3]-T[j+1])/2) * A * (T_e-(T[j+3]-T[j+1])/2.0)
        q.append(q_a1)
        
        # Even Region (Heat input)
        q_w2 = m_w * Cp('w',(T[j+4]-T[j+2])/2) * (T[j+4]-T[j+2]) + q_receiver
        q.append(q_w2)
        q_a2 = m_a * Cp('a',(T[j+5]-T[j+3])/2) * (T[j+5]-T[j+3]) + q_module
        q.append(q_a2)
        j = j + 4
    return (q)
        
T = GenTemperatureArray(m)
print T
q = Residue(m,T)
print q