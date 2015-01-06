# Written by Justin Shultz from Rensselaer Polytechnic Institute 
# for Assad Oberai, with his assistance too

import numpy
import math

def Cp(liquid, T):
    # This function calculates the Specific Heat of Water and Air based on Temperature
    # The equation for Water was a 6 degree polynomial calculated using a temperature range from 5-105C
    # The specifc heat of air does not change over the temperature range of the system and can be set constant
    if liquid == "w":
        # Constants determined from a 5 degree fit polynomial using MatLab CurveFit Toolbox
        # R = 0.9986 & SSE = 5.556e-06
        p1 =   -0.001201
        p2 =     0.00218
        p3 =    0.001323 
        p4 =     0.00634
        p5 =     0.01152
        p6 =       4.183
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
    
    if interface == 'p':
        k = k('w', T)
        L = 3*10**-3
        Nu = 3.66 # for Reynold numbers < 2300 (laminar) and constant wall temperature
    elif interface == 'pa':
        k = k('a', T)
        L = 0.3 # NOT LENGTH OF PIPE 
        Re  = rho_a(T) * m_a * radius_OD / viscosity_a(T)
        Nu = 0.037*Re**(4.0/5.0)*Pr_a(T)**(1.0/3.0)
    elif interface == 'i':
        k = k('a', T)
        L = 3.0
        Re  = rho_a(T) * m_a * radius_cavity / viscosity_a(T)
        Nu = 0.037*Re**(4.0/5.0)*Pr_a(T)**(1.0/3.0)
    elif interface == 'e':
        k = k('a', T)
        L = 3.0
        Re  = rho_a(T) * m_a * radius_cavity / viscosity_a(T)
        Nu = 0.037*Re**(4.0/5.0)*Pr_a(T)**(1.0/3.0)
    h = k / L * Nu # L is the diameter or boundary layer thickness
    return h

def R(interface,T):
    if interface == "wa":  #water to air
        # Convenction_interior + Silicone Pipe + Insulation + Convection_exterior
        # A is cross sectional area
        R_conv_wt = 1/(h('p',T) * math.pi * (3*10**-3)**2)    # Convection from water to tubing
        R_cond_t = math.pi*((3+1.5)-3)*10**-3 / (0.145 * math.pi * ((3+1.5)*10**-3)**2)    # Conduction through tubing
        R_cond_ins = math.pi*((19.05)-4.5)*10**-3/(0.037 * math.pi * ((ins_pipe+3+1.5)*10**-3)**2)   # Conduction through insulation
        R_conv_ta = 1/(h('pa',T) * math.pi * ((3+1.5+ins_pipe)*10**-3)**2)    # Convenction from tubing to air
        Resistance = R_conv_wt + R_cond_t + R_cond_ins + R_conv_ta
    
    if interface == 'i': # Resistance of the interior wall
        # L/k = thickness/conductivity
        # Double Layer Window with Argon Gap
        FirstLayer = (6*10**-3)/(1.05) # First layer = first layer of glass
        SecondLayer = (6*10**-3)/(0.016) # Second layer = argon gap
        ThirdLayer = (6*10**-3)/(1.05) # Third layer = second layer of glass
        Resistance = FirstLayer + SecondLayer + ThirdLayer
    
    if interface == 'e': # Resistance of the exterior wall
        # L/k = thickness/conductivity
        # Single layer exterior optically transparent glass layer for ideal optical performance on ICSolar
        FirstLayer = (6*10**-3)/(1.05) # First layer = first layer of glass
        Resistance = FirstLayer
    return Resistance

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
        q_w1 = (m_w * Cp('w',(T[j+2]-T[j])/2) * (T[j+2]-T[j])) - (As_pipe * ((T[j+3]-T[j+1])/2.0-(T[j+2]-T[j])/2.0) / R('wa',(T[j+3]-T[j+1])/2.0-(T[j+2]-T[j])/2))
        q[j] = q_w1
        # Calculate air temperature for region 1 (Transfer with water, interior and exterior)
        q_a1 = m_a * Cp('a',(T[j+3]-T[j+1])/2) * (T[j+3]-T[j+1]) + (As_pipe * ((T[j+3]-T[j+1])/2.0 - (T[j+2]-T[j])/2.0)/R('wa',(T[j+3]-T[j+1])/2.0-(T[j+3]-T[j+1])/2.0)) - h('i',T_i-(T[j+3]-T[j+1])/2) * A_wall * (T_i-(T[j+3]-T[j+1])/2.0) - h('e',T_e-(T[j+3]-T[j+1])/2) * A_wall * (T_e-(T[j+3]-T[j+1])/2.0)
        q[j+1] = q_a1
        
        # Even Region (Heat input)
        q_w2 = m_w * Cp('w',(T[j+4]-T[j+2])/2) * (T[j+4]-T[j+2]) + q_receiver
        q[j+2] = q_w2
        q_a2 = m_a * Cp('a',(T[j+5]-T[j+3])/2) * (T[j+5]-T[j+3]) + q_module
        q[j+3] = q_a2
        j = j + 4
    return (q)
   
m = 3 # Number of Modules
As_pipe = 1.0 # Area m^2
A_wall = 3*2 # Area of glass
ins_pipe = 19.05 # Thickness of the insulation around the silicone pipes

q_receiver = 200.0 # Heat flow into water from Module Heat Receiver
q_module = 50.0 # Heat flow into air from Heat Loss from the Module

m_w = 4.0 *10**(-6) * rho_w((13+90)/2) # Mass flowrate of water
m_a = 4.0 # Flowrate of air

T_i = 20.0 # Temperature on the interior of the building, degrees C
T_e = 25.0 # Temperature on the exterior of the building, degrees C
T_ai = 20.0 # Inlet temperature of air, degrees C
T_wi = 13.0 # Inlet temperature of water, degrees C     
     
T = GenTemperatureArray(m)
print T
q = Residue(m,T)
print q