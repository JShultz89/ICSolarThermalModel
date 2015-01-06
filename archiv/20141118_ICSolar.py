# Written by Justin Shultz from Rensselaer Polytechnic Institute 
# for Assad Oberai, with his assistance 

import numpy
import math

def Cp(liquid, T):
    # This function calculates the Specific Heat of Water and Air based on Temperature
    # The equation for Water was a 6 degree polynomial calculated using a temperature range from 5-105C
    # The specifc heat of air does not change over the temperature range of the system and can be set constant
    if liquid == "w":
        # Constants determined from a 5 degree fit polynomial using MatLab CurveFit Toolbox
        # R = 0.9986 & SSE = 5.556e-06
        p1 =  -4.178e-11
        p2 =   1.384e-08
        p3 =  -1.737e-06
        p4 =   0.0001115
        p5 =   -0.003429
        p6 =       4.218
        Cp = p1*T**5 + p2*T**4 + p3*T**3 + p4*T**2 + p5*T + p6
    elif liquid == "a":
        # The specifc heat of air does not change over the temperature region of the system
        Cp = 1.005
    return Cp # Return the specific heat

def rho_a(T):
    p1 =    1.75e-05
    p2 =    -0.00483 
    p3 =       1.293 
    rho = p1*T**2 + p2*T + p3
    return rho
    
def rho_w(T):
    p1 =   -0.003416
    p2 =    -0.09298
    p3 =        1001
    rho = p1*T**2 + p2*T + p3
    return rho

def h(interface, T): 
    # TODO: Figure out characteristic length
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
            # Constants determined from a 1 degree fit polynomial using MatLab CurveFit Toolbox
            # R = 1
            p1 =       7e-05
            p2 =      0.0243
            k = p1*T + p2
        return k # Return the specific heat
        
    def viscosity_a(T):
        p1 =     7.5e-05
        p2 =      0.0888
        p3 =        13.3
        viscosity = p1*T**2 + p2*T + p3
        return viscosity 
        
    def Pr_a(T):
        p1 =  -4.705e-19
        p2 =     -0.0001
        p3 =       0.715 
        Pr = p1*T**2 + p2*T + p3
        return Pr
    
    radius_OD = 3./4. * 0.0254 # Outside Diameter Radius of Pipe with Insulation
    radius_cavity = 4 * 0.3*0.3 / (0.3*4) / 2  # Hydraulic diameter assuming cavity is 0.3 meters by 0.3 m (cross section)
    
    if mode == 'forced':
        # Characteristic length discussion: http://physics.stackexchange.com/questions/44185/what-is-the-characteristic-length-of-a-cylinder
        if interface == 'p':
            k = k('w', T)
            L = 0.3 # Diameter
            Nu = 3.66 # Nu_D for Reynold numbers < 2300 (laminar) and constant wall temperature
        elif interface == 'pa':
            k = k('a', T)
            L = 0.3 * m # Length of system
            Re  = rho_a(T) * v_a * 0.3*m / viscosity_a(T) # Re = velocity*Desnity*Length/viscosity
            Nu = 0.037*Re**(4.0/5.0)*Pr_a(T)**(1.0/3.0) 
        elif interface == 'i':
            k = k('a', T)
            L = 0.3*m
            Re = rho_a(T) * v_a * 0.3*m / viscosity_a(T)
            Nu = 0.037*Re**(4.0/5.0)*Pr_a(T)**(1.0/3.0)
        elif interface == 'e':
            k = k('a', T)
            L = 0.3*m
            Re  = rho_a(T) * v_a * 0.3*m / viscosity_a(T)
            Nu = 0.037*Re**(4.0/5.0)*Pr_a(T)**(1.0/3.0)
        h = k / L * Nu # L the characteristic length that should match the Nu L term
    return h

def R(interface, T):
    if interface == "wa":  #water to air
        # Resistance = Convenction_interior + Silicone Pipe + Insulation + Convection_exterior
        # A is cross sectional area
        R_conv_wt = 1.0/(2* math.pi * (3*10**-3) * 0.3 * h('p',(90+13)/2))    # Convection from water to tubing
        R_cond_t = ((3-1.5)*10**-3) / (2 * math.pi * (3+1.5)*10**-3 * 0.3 * 0.145)    # Conduction through tubing
        R_cond_ins = ((ins_pipe)*10**-3) / (2.0 * math.pi * (3+1.5+ins_pipe)*10**-3 * 0.3 * 0.037)    # Conduction through insulation
        R_conv_ta = 1.0/(2 * math.pi * (3+1.5+ins_pipe)*10**-3 * 0.3 * h('pa',(30+13)/2))    # Convenction from tubing to air
        Resistance = R_conv_wt + R_cond_t + R_cond_ins + R_conv_ta
    
    if interface == 'i': # Resistance of the interior wall
        # L/kA = thickness/(conductivity*Area)
        # Double Layer Window with Argon Gap
        # Resistance = First Layer of Glass + Argon Gap + Second Layer of Glass
        FirstLayer = (6*10**-3)/(1.05 * A_wall) # First layer = first layer of glass
        SecondLayer = (6*10**-3)/(0.016 * A_wall) # Second layer = argon gap
        ThirdLayer = (6*10**-3)/(1.05 * A_wall) # Third layer = second layer of glass
        Resistance = FirstLayer + SecondLayer + ThirdLayer
    
    if interface == 'e': # Resistance of the exterior wall
        # L/kA = thickness/(conductivity*area)
        # Single layer exterior optically transparent glass layer for ideal optical performance on ICSolar
        FirstLayer = (6*10**-3)/(1.05 * A_wall) # First layer = first layer of glass
        Resistance = FirstLayer
    return Resistance

def GenTemperatureArray(m):
    # This function creates an array of temperature values (Even is air, Odd is water)
    # The values are an initial guess for use in the residue function
    count = 0
    T = numpy.empty(m*4+2, dtype=float)
    T_w_tmp = T_wi # Temperary value for the water 
    T_a_tmp = T_ai # Temperary value for the water
    guessIncrement = 5 # Increases temperature in each stage by number, may not be necessary
    while count < 4*m+2:
        if count%2 == 0:
            T_tmp = T_w_tmp
            T[count] = T_tmp # Add the first air temperature to the array
            T_w_tmp += guessIncrement # Increase each air temperature by 5 degrees 
        else:
            T_tmp = T_a_tmp
            T[count] = T_tmp # Add the first water temeperature to the end of the array
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
        q_w1 = m_w * Cp('w', ((T[j+2]-T[j])/2)*(T[j+2]-T[j])) - ( ((T[j+3]-T[j+1])/2.0-(T[j+2]-T[j])/2.0) / R('wa',(T[j+3]-T[j+1])/2.0-(T[j+2]-T[j])/2) )
        q[j] = q_w1
        # Calculate air temperature for region 1 (Transfer with water, interior and exterior)
        q_a1 = m_a * Cp('a',(T[j+3]-T[j+1])/2) * (T[j+3]-T[j+1]) + (((T[j+3]-T[j+1])/2.0 - (T[j+2]-T[j])/2.0)/R('wa',(T[j+3]-T[j+1])/2.0-(T[j+2]-T[j])/2.0)) - (T_int-(T[j+3]-T[j+1])/2.0)/R('i',T_int-(T[j+3]-T[j+1])/2) - (T_ext-(T[j+3]-T[j+1])/2.0)/R('e',T_ext-(T[j+3]-T[j+1])/2)
        q[j+1] = q_a1
        
        # Even Region (Heat input)
        q_w2 = m_w * Cp('w',(T[j+4]-T[j+2])/2) * (T[j+4]-T[j+2]) + q_receiver
        q[j+2] = q_w2
        q_a2 = m_a * Cp('a',(T[j+5]-T[j+3])/2) * (T[j+5]-T[j+3]) + q_module
        q[j+3] = q_a2
        j = j + 4
    return (q)

def Differentiation(q,T):
    dRdT = numpy.zeros((m*4,m*4), dtype=float)
    j = 0
    for i in range(0,m):
        ##### Region 1: Heat Balance in the Water #####
        # Derivative of R1 (water balance) with respect to T[j] (T[0])
        dR1dT0_w = - m_w * Cp('w', ((T[j+2]-T[j])/2)) - 1.0/2.0 / R('wa',(T[j+3]-T[j+1])/2.0-(T[j+2]-T[j])/2.0)
        # Outside of Matrix
        
        # Derivative of R1 (water balance) with respect to T[j+1] (T[1]) 
        dR1dT1_w = 1.0/2.0 / R('wa',(T[j+3]-T[j+1])/2.0-(T[j+2]-T[j])/2)
        # Outside of Matrix
    
        # Derivative of R1 (water balance) with respect to T[j+2] (T[2])
        dR1dT2_w = m_w * Cp('w', ((T[j+2]-T[j])/2)) + 1.0/2.0 / R('wa',(T[j+3]-T[j+1])/2.0-(T[j+2]-T[j])/2.0)
        dRdT[j,j] = T[j+2] - dR1dT2_w/q[j]
    
        # Derivative of R1 (water balance) with respect to T[j+3] (T[3]) 
        dR1dT3_w = - 1.0/2.0 / R('wa',(T[j+3]-T[j+1])/2.0-(T[j+2]-T[j])/2)
        dRdT[j,j+1] = T[j+3] - dR1dT3_w/q[j+1]
    
        ##### Region 1: Heat Balance in the Air #####
        # Derivative of R1 (air balance) with respect to T[j] (T[0])
        dR1dT0_a = (-(-1.0)/2.0)/R('wa',(T[j+3]-T[j+1])/2.0-(T[j+2]-T[j])/2.0)
        # Outside of Matrix
    
        # Derivative of R1 (air balance) with respect to T[j+1] (T[1])
        dR1dT1_a = - m_a * Cp('a',(T[j+3]-T[j+1])/2) + (((-1.0)/2.0)/R('wa',(T[j+3]-T[j+1])/2.0 - (T[j+2]-T[j])/2.0)) - (-(-1.0)/2.0)/R('i',T_int-(T[j+3]-T[j+1])/2.0) - (-(-1.0)/2.0)/R('e',T_ext-(T[j+3]-T[j+1])/2)
        # Outside of Matrix    
        
        # Derivative of R1 (air balance) with respect to T[j+2] (T[2])
        dR1dT2_a = ((1)/2.0)/R('wa',(T[j+3]-T[j+1])/2.0-(T[j+2]-T[j])/2.0)
        dRdT[j+1,j] = T[j+2] - dR1dT2_a/q[j]
    
        # Derivative of R1 (air balance) with respect to T[j+3] (T[3])
        dR1dT3_a = m_a * Cp('a',(T[j+3]-T[j+1])/2) + 1.0/2.0/R('wa',(T[j+3]-T[j+1])/2.0-(T[j+2]-T[j])/2.0) - (-(1.0)/2.0)/R('i',T_int-(T[j+3]-T[j+1])/2) - (-(1.0)/2.0)/R('e',T_ext-(T[j+3]-T[j+1])/2.0)
        dRdT[j+1,j+1] = T[j+3] - dR1dT3_a/q[j+1]
        
        
        
        ##### Region 2: Heat Balance in the Water #####
        # Derivative of R2 (water balance) with respect to T[j] (T[0])
        dR2dT0_w = 0
        # Outside of Matrix  
        
        # Derivative of R2 (water balance) with respect to T[j+1] (T[1]) 
        dR2dT1_w = 0
        # Outside of Matrix
    
        # Derivative of R2 (water balance) with respect to T[j+2] (T[2])
        dR2dT2_w = - m_w * Cp('w',(T[j+4]-T[j+2])/2)
        dRdT[j+2,j] = T[j+2] - dR2dT2_w/q[j]
    
        # Derivative of R2 (water balance) with respect to T[j+3] (T[3]) 
        dR2dT3_w = 0
    
        # Derivative of R2 (water balance) with respect to T[j+4] (T[4]) 
        dR2dT4_w = m_w * Cp('w',(T[j+4]-T[j+2])/2)
        dRdT[j+2,j+3] = T[j+4] - dR2dT4_w/q[j+2]
        
        # Derivative of R2 (water balance) with respect to T[j+5] (T[5]) 
        dR2dT5_w = 0
        
        ##### Region 2: Heat Balance in the Air #####
        # Derivative of R2 (air balance) with respect to T[j] (T[0])
        dR2dT0_a = 0
        # Outside of Matrix  
        
        # Derivative of R2 (air balance) with respect to T[j+1] (T[1]) 
        dR2dT1_a = 0
        # Outside of Matrix
    
        # Derivative of R2 (air balance) with respect to T[j+2] (T[2])
        dR2dT2_a = 0
    
        # Derivative of R2 (air balance) with respect to T[j+3] (T[3]) 
        dR2dT3_a = - m_a * Cp('a',(T[j+5]-T[j+3])/2)
        dRdT[j+3,j+1] = T[j+3] - dR2dT3_a/q[j+2]
    
        # Derivative of R2 (air balance) with respect to T[j+4] (T[4]) 
        dR2dT4_a = 0
        
        # Derivative of R2 (air balance) with respect to T[j+5] (T[5]) 
        dR2dT5_a = m_a * Cp('a',(T[j+5]-T[j+3])/2)
        dRdT[j+3,j+3] = T[j+5] - dR2dT5_a/q[j+3]
        
        j = j + 4
    
    return dRdT


m = 1 # Number of Modules
ins_pipe = 19.05 # Thickness of the insulation around the silicone pipes
As_pipe = 2 * math.pi * ((3+1.5+ins_pipe)*10**-3) * 0.3 # Area m^2
Height_wall = 3
Width_wall = 0.3
A_wall = Height_wall * Width_wall # Area of glass
P_wall = 2*Height_wall + 2*Width_wall # Perimeter of the wall

q_receiver = 8.0 # Heat flow into water from Module Heat Receiver
q_module = 3.0 # Heat flow into air from Heat Loss from the Module

mode = 'forced' # Defines whether the air flow is 'forced' or 'natural'
m_w = 4.0 *10**(-6) * rho_w((13+90)/2) # Mass flowrate of water = VolumetricFlowrate * DensityWater
m_a = 4.0 * rho_a((13+30)/2) # Mass Flowrate of air = VolumetricFlowrate * DensityAir
v_a = m_a / (0.4*0.4 * rho_a((30-13)/2)) # Flow velocity of air = mass flow rate / (rho of air * cross sectional area)

T_int = 20.0 # Temperature on the interior of the building, degrees C
T_ext = 25.0 # Temperature on the exterior of the building, degrees C
T_ai = 20.0 # Inlet temperature of air, degrees C
T_wi = 13.0 # Inlet temperature of water, degrees C     
     
T = GenTemperatureArray(m)
print T
q = Residue(m,T)
print q
T_new = Differentiation(q,T)
print T_new