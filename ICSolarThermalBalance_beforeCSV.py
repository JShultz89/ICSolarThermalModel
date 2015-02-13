# Written by Justin Shultz from Rensselaer Polytechnic Institute 
# for Assad Oberai, with his assistance 

import numpy
import math
import time

numpy.set_printoptions(linewidth=150)

def Cp(liquid, T):
    # This function calculates the fSpecific Heat of Water and Air based on Temperature
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
    print 'The specific heat of ' ,str(liquid), " is ", Cp
    return Cp # Return the specific heat

def rho_a(T): # kg/m^3
    p1 =    1.75e-05
    p2 =    -0.00483 
    p3 =       1.293 
    rho = p1*T**2 + p2*T + p3
    print "Density of ", str(rho_a), " is ", rho
    return rho
    
def rho_w(T):
    p1 =   -0.003416
    p2 =    -0.09298
    p3 =        1001
    rho = p1*T**2 + p2*T + p3
    print "Density of", str(rho_w), " is ", rho
    return rho

def h(interface, T): 
    # This function returns the convection heat transfer coefficient for specific temperatures
    # Interface, at which exchange is occuring. Options: w (water interior of pipe), pa (pipe to air), i (interior), e (exterior)
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
        viscosity = (p1*T**2 + p2*T + p3)*10**(-6)
        return viscosity 
        
    def Pr_a(T):
        p1 =  -4.705e-19
        p2 =     -0.0001
        p3 =       0.715 
        Pr = p1*T**2 + p2*T + p3
        return Pr

    if mode == 'forced':
        # Characteristic length discussion: http://physics.stackexchange.com/questions/44185/what-is-the-characteristic-length-of-a-cylinder
        if interface == 'w':
            k = k('w', T)
            L = diameter_inside # Diameter *NOTE* Check this number again had 0.3 before, that doesn't make sense because it should be milimeters
            Nu = 3.66 # Nu_D for Reynold numbers < 2300 (laminar) and constant wall temperature. Could use Nu_D = 4.36 for constant heat transfer.
        elif interface == 'pa':
            k = k('a', T)
            L = tubing_length * m # Length of system
            Re  = v_a * L / viscosity_a(T) # Re = velocity*Desnity*Length/viscosity
            Nu = 0.037*Re**(4.0/5.0)*Pr_a(T)**(1.0/3.0) 
        elif interface == 'i':
            k = k('a', T)
            L = tubing_length * m
            Re = v_a * L / viscosity_a(T)
            Nu = 0.0296*Re**(4.0/5.0)*Pr_a(T)**(1.0/3.0)
        elif interface == 'e':
            k = k('a', T)
            L = tubing_length * m
            Re  = v_a * L / viscosity_a(T)
            Nu = 0.0296*Re**(4.0/5.0)*Pr_a(T)**(1.0/3.0)
        h = k / L * Nu # L the characteristic length that should match the Nu L term
    print "Convection heat transfer coefficient of ", str(interface), " is ", h
    return h

def R(interface):
    if interface == "pipe":  # water to air
        # Resistance = Convenction_interior + Silicone Pipe + Insulation + Convection_exterior
        # A is cross sectional area
        R.R_conv_wt = 1.0/(2.0 * math.pi * diameter_inside * tubing_length * h('w', (T[0]+T[m*4])/2.0 ))    # Convection from water to tubing
        R.R_cond_t = math.log((diameter_tubing)/(diameter_inside)) / (2 * math.pi * tubing_length * cond_tubing)    # Conduction through tubing
        R.R_cond_ins = math.log((diameter_insulation)/(diameter_tubing)) / (2.0 * math.pi * tubing_length * cond_insulation)    # Conduction through insulation
        R.R_conv_ta = 1.0/(2 * math.pi * diameter_insulation * tubing_length * h('pa',(T[0+1]+T[m*4+1])/2+(T[0+1]+T[m*4+1])/2))    # Convenction from tubing to air
        Resistance = R.R_conv_wt + R.R_cond_t + R.R_cond_ins + R.R_conv_ta
    
    if interface == 'int': # Resistance of the interior wall
        # L/kA = thickness/(conductivity*Area)
        # Double Layer Window with Argon Gap
        # Resistance = First Layer of Glass + Argon Gap + Second Layer of Glass
        Interior_Convection = 1/(A_wall * h('i',(((T[0+1]+T[m*4+1])/2)+T_int)/2))
        FirstLayer = (6*10**(-3))/(1.05 * A_wall) # First layer = first layer of glass
        SecondLayer = (6*10**-3)/(0.016 * A_wall) # Second layer = argon gap
        ThirdLayer = (6*10**-3)/(1.05 * A_wall) # Third layer = second layer of glass
        Resistance = FirstLayer + SecondLayer + ThirdLayer + Interior_Convection
    
    if interface == 'ext': # Resistance of the exterior wall
        # L/kA = thickness/(conductivity*area)
        # Single layer exterior optically transparent glass layer for ideal optical performance on ICSolar
        Interior_Convection = 1/(A_wall * h('e',(((T[0+1]+T[m*4+1])/2)+T_ext)/2))
        FirstLayer = (6*10**-3)/(1.05 * A_wall) # First layer = first layer of glass
        Resistance = FirstLayer + Interior_Convection
    print "The resistance of ", str(interface), " is ", Resistance
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
            T_w_tmp += 2 # Increase each air temperature by 5 degrees 
        else:
            T_tmp = T_a_tmp
            T[count] = T_tmp # Add the first water temeperature to the end of the array
            T_a_tmp += 5 # Increase each water temperature by 5 degrees
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
        q_w1 = m_w * Cp('w', ((T[j+2]+T[j])/2)*(T[j+2]-T[j])) - ((T[j+3]-T[j+2]) / R('pipe'))
        q[j] = q_w1
        # Calculate air temperature for region 1 (Transfer with water, interior and exterior)
        q_a1 = m_a * Cp('a',(T[j+3]+T[j+1])/2) * (T[j+3]-T[j+1]) + ((T[j+3]-T[j+2]) / R('pipe')) - (T_int-T[j+3])/R('int') - (T_ext-T[j+3])/R('ext')
        q[j+1] = q_a1
        
        # Even Region (Heat input)
        q_w2 = m_w * Cp('w',(T[j+4]+T[j+2])/2) * (T[j+4]-T[j+2]) - q_receiver
        q[j+2] = q_w2
        q_a2 = m_a * Cp('a',(T[j+5]+T[j+3])/2) * (T[j+5]-T[j+3]) - q_module
        q[j+3] = q_a2
        j = j + 4
    return (q)

def Differentiation():
    dRdT = numpy.zeros((m*4,m*4), dtype=float)
    j = 0
    for i in range(0,m):
            ##### Region 1: Heat Balance in the Water #####
            # q_w1 = m_w * Cp('w', ((T[j+2]+T[j])/2)*(T[j+2]-T[j])) - ((T[j+3]-T[j+2]) / R('pipe')) 
        # Derivative of R1 (water balance) with respect to T[j] (T[0])
        dR1dT0_w = - m_w * Cp('w', (T[j+2]+T[j])/2)
        if j == 0:
            #Do nothing
            pass
        else:
            dRdT[j,j-2] = dR1dT0_w
        
        # Derivative of R1 (water balance) with respect to T[j+1] (T[1]) 
        dR1dT1_w = 0
        if j == 0:
            #Do nothing
            pass
        else:
            dRdT[j,j-1] = dR1dT1_w
    
        # Derivative of R1 (water balance) with respect to T[j+2] (T[2])
        dR1dT2_w = m_w * Cp('w', ((T[j+2]+T[j])/2)) + 1.0/R('pipe')
        dRdT[j,j] = dR1dT2_w
    
        # Derivative of R1 (water balance) with respect to T[j+3] (T[3]) 
        dR1dT3_w = - 1.0/R('pipe')
        dRdT[j,j+1] = dR1dT3_w
    
    
    
    
            ##### Region 1: Heat Balance in the Air #####
            # q_a1 = m_a * Cp('a',(T[j+3]+T[j+1])/2) * (T[j+3]-T[j+1]) + ((T[j+3]-T[j+2]) / R('pipe')) - (T_int-T[j+3])/R('int') - (T_ext-T[j+3])/R('ext')
        # Derivative of R1 (air balance) with respect to T[j] (T[0])
        dR1dT0_a = 0
        if j == 0:
            #Do nothing
            pass
        else:
            dRdT[j+1,j-2] = dR1dT0_a
    
        # Derivative of R1 (air balance) with respect to T[j+1] (T[1])
        dR1dT1_a = - m_a * Cp('a',(T[j+3]+T[j+1])/2) 
        if j == 0:
            #Do nothing
            pass
        else:
            dRdT[j+1,j-1] = dR1dT1_a    
        
        # Derivative of R1 (air balance) with respect to T[j+2] (T[2])
        dR1dT2_a = -1.0/R('pipe')
        dRdT[j+1,j] = dR1dT2_a
    
        # Derivative of R1 (air balance) with respect to T[j+3] (T[3])
        dR1dT3_a = m_a * Cp('a',(T[j+3]+T[j+1])/2.0) + 1.0/R('pipe') + 1.0/R('int') + 1.0/R('ext')
        dRdT[j+1,j+1] = dR1dT3_a
        
        
        
            ##### Region 2: Heat Balance in the Water #####
            # q_w2 = m_w * Cp('w',(T[j+4]+T[j+2])/2) * (T[j+4]-T[j+2]) + q_receiver
        # Derivative of R2 (water balance) with respect to T[j] (T[0])
        dR2dT0_w = 0
        
        # Derivative of R2 (water balance) with respect to T[j+1] (T[1]) 
        dR2dT1_w = 0
    
        # Derivative of R2 (water balance) with respect to T[j+2] (T[2])
        dR2dT2_w = - m_w * Cp('w',(T[j+4]+T[j+2])/2)
        dRdT[j+2,j] = dR2dT2_w
    
        # Derivative of R2 (water balance) with respect to T[j+3] (T[3]) 
        dR2dT3_w = 0
    
        # Derivative of R2 (water balance) with respect to T[j+4] (T[4]) 
        dR2dT4_w = m_w * Cp('w',(T[j+4]+T[j+2])/2)
        dRdT[j+2,j+2] = dR2dT4_w
        
        # Derivative of R2 (water balance) with respect to T[j+5] (T[5]) 
        dR2dT5_w = 0
        
        
            ##### Region 2: Heat Balance in the Air #####
            # q_a2 = m_a * Cp('a',(T[j+5]+T[j+3])/2) * (T[j+5]-T[j+3]) + q_module
        # Derivative of R2 (air balance) with respect to T[j] (T[0])
        dR2dT0_a = 0
        # Outside of Matrix  
        
        # Derivative of R2 (air balance) with respect to T[j+1] (T[1]) 
        dR2dT1_a = 0
        # Outside of Matrix
    
        # Derivative of R2 (air balance) with respect to T[j+2] (T[2])
        dR2dT2_a = 0
    
        # Derivative of R2 (air balance) with respect to T[j+3] (T[3]) 
        dR2dT3_a = - m_a * Cp('a',(T[j+5]+T[j+3])/2)
        dRdT[j+3,j+1] = dR2dT3_a
    
        # Derivative of R2 (air balance) with respect to T[j+4] (T[4]) 
        dR2dT4_a = 0
        
        # Derivative of R2 (air balance) with respect to T[j+5] (T[5]) 
        dR2dT5_a = m_a * Cp('a',(T[j+5]+T[j+3])/2)
        dRdT[j+3,j+3] = dR2dT5_a
        
        j = j + 4
    
    return dRdT


def UpdatedTemp(T):    
    def Increment(dQ, q):
        # Using Newton Method
        # Increment = q/dQ = dQ_inv * q
        # T_new = T_old - dQ_inv * q = T_old - q/dQ
        q_matrix = numpy.matrix(q)
        dQ_inv = numpy.linalg.inv(dQ)
        Increment = dQ_inv * numpy.transpose(q_matrix)
        Increment = numpy.squeeze(numpy.asarray(Increment))
        return Increment
        
    print "Inverse Derivative Matrix"
    print numpy.linalg.inv(dQ)
    print 
    Inc = Increment(dQ, q)
    print "This is the initial increment"
    print Inc
    print 
    
    T_new = T[2:]
    
    count = 1
    
    # while count < 2: #abs(numpy.amax(Inc)) > 10**(-4):
    T_new = T_new - Inc 
    T_local = numpy.concatenate([T[:2], T_new])
    #print "This is the temperature after %d increment" %(count)
    #print T_local
    #print
    
    #q_new = Residue(m,T_local)
    #print "Updating heat balance with new temperature"
    #print q_new
    #print
    
    #Inc = Increment(dQ, q_new)
    #print "Updated increment"
    #print Inc
    #print
    count += 1
    
    T_final = numpy.concatenate([T[:2], T_new])
    return T_final
    

m = 1 # Number of Modules

T_int = 22.5 # Temperature on the interior of the building, degrees [C]
T_ext = 25.0 # Temperature on the exterior of the building, degrees [C]
T_ai = 20.0 # Inlet temperature of air, degrees [C]
T_wi = 13.0 # Inlet temperature of water, degrees [C]  

T = GenTemperatureArray(m)
print "Initial Temperature Array Guess"
print T
print 

diameter_inside = 3.0*10**(-3)
diameter_tubing = diameter_inside + 1.675*10**(-3)
ins_thickness = 9.525*10**(-3) # Thickness of the insulation around the silicone pipes
diameter_insulation = diameter_tubing + ins_thickness
As_pipe = 2 * math.pi * (diameter_insulation) * 0.3 # Area m^2
cond_tubing = 0.145 # Conduction value of silicon tubing
cond_insulation = 0.037 # Conduction value of silicon insulation

Height_wall = 0.3
Width_wall = 0.3
Depth_caviety = 0.3
tubing_length = 0.3
A_wall = Height_wall * Width_wall # Area of glass
P_wall = 2*Height_wall + 2*Width_wall # Perimeter of the wall
radius_cavity = 4 * 0.3*0.3 / (0.3*4) / 2  # Hydraulic diameter assuming cavity is 0.3 meters by 0.3 m (cross section)


q_receiver = 0 # 8.0*10**(-3) # Heat flow into water from Module Heat Receiver
q_module = 0 # 3.0*10**(-3) # Heat flow into air from Heat Loss from the Module

mode = 'forced' # Defines whether the air flow is 'forced' or 'natural'
m_w = 8.5*10**(-7) * rho_w((T[0]+T[m*4])/2) # Mass flowrate of water = VolumetricFlowrate * DensityWater
v_a = 2.0 # Flow velocity [m/s]
m_a = v_a * rho_a((13+30)/2) * Depth_caviety*Width_wall # Mass Flowrate of air [kg/s] = velocity [m/s] * DensityAir [kg/m^3] * cross section [m^2]   
     

q = Residue(m,T)
print "Residual heat balance array:"
print q
print 
dQ = Differentiation()
print "Derivative Matrix:"
print dQ
print
T_final = UpdatedTemp(T)
print "Final Temperatures"
print T_final