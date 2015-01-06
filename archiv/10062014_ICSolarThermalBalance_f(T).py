# Written by Justin Shultz from Rensselaer Polytechnic Institute 
# for Assad Oberai, with his assistance too

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
T_wi = 13.0 # Inlet temperature of water, degrees C
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

# Finding solution in which i (for water) and i+1 (for air) returns a 0 result
#def residue
def Residue(m,T):
    q = []
    j = 0
    for i in range(0,m):
        print "Running module %d" %(i+1)
        # Odd region (transfer between water, air, interior, and exterior)
        # Calculate water temperature for region 1 (Transfer with air)
        q_w1 = m_w * C_pw * (T[j+2]-T[j]) - h_wa * A * ((T[j+3]-T[j+1]/2.0-(T[j+2]-T[j])/2.0))
        q.append(q_w1)
        # Calculate air temperature for region 1 (Transfer with water, interior and exterior)
        q_a1 = m_a * C_pa * (T[j+3]-T[j+1]) + h_wa * A * ((T[j+3]-T[j+1])/2.0-(T[j+2]-T[j])/2.0)-h_i*A*(T_i-(T[j+3]-T[j+1])/2.0)-h_e*A*(T_e-(T[j+3]-T[j+1])/2.0)
        q.append(q_a1)
        
        # Even Region (Heat input)
        q_w2 = m_w * C_pw * (T[j+4]-T[j+2]) + q_receiver
        q.append(q_w2)
        q_a2 = m_a * C_pa * (T[j+5]-T[j+3]) + q_module
        q.append(q_a2)
        j = j + 4
    return (q)
        
T = GenTemperatureArray(3)
print T
q = Residue(m,T)
print q