def Differentiation(q,T):
   # Calculate the differentiation of the heat equations in the residue
   
   ##### Region 1: Heat Balance in the Water #####
   # Derivative of R1 (water balance) with respect to T[j] (T[0])
   dR1dT0_w = m_w * Cp('w', ((T[j+2]-T[j])/2)*(T[j+2]-0)) - ( ((T[j+3]-T[j+1])/2.0-(T[j+2]-0)/2.0) / R('wa',(T[j+3]-T[j+1])/2.0-(T[j+2]-T[j])/2) )

   # Derivative of R1 (water balance) with respect to T[j+1] (T[1])
   dR1dT1_w = m_w * Cp('w', ((T[j+2]-T[j])/2)*(T[j+2]-T[j])) - ( ((T[j+3]-0)/2.0-(T[j+2]-T[j])/2.0) / R('wa',(T[j+3]-T[j+1])/2.0-(T[j+2]-T[j])/2) )

   # Derivative of R1 (water balance) with respect to T[j+2] (T[2])
   dR1dT2_w = m_w * Cp('w', ((T[j+2]-T[j])/2)*(0-T[j])) - ( ((T[j+3]-T[j+1])/2.0-(0-T[j])/2.0) / R('wa',(T[j+3]-T[j+1])/2.0-(T[j+2]-T[j])/2) )

   # Derivative of R1 (water balance) with respect to T[j+3] (T[3])
   dR1dT3_w = m_w * Cp('w', ((T[j+2]-T[j])/2)*(T[j+2]-T[j])) - ( ((0-T[j+1])/2.0-(T[j+2]-T[j])/2.0) / R('wa',(T[j+3]-T[j+1])/2.0-(T[j+2]-T[j])/2) )

   ##### Region 1: Heat Balance in the Air #####
   # Derivative of R1 (air balance) with respect to T[j] (T[0])
   dR1dT0_a = m_a * Cp('a',(T[j+3]-T[j+1])/2) * (T[j+3]-T[j+1]) + (((T[j+3]-T[j+1])/2.0 - (T[j+2]-0)/2.0)/R('wa',(T[j+3]-T[j+1])/2.0-(T[j+2]-T[j])/2.0)) - (T_int-(T[j+3]-T[j+1])/2.0)/R('i',T_int-(T[j+3]-T[j+1])/2) - (T_ext-(T[j+3]-T[j+1])/2.0)/R('e',T_ext-(T[j+3]-T[j+1])/2)

   # Derivative of R1 (air balance) with respect to T[j+1] (T[1])
   dR1dT1_a = m_a * Cp('a',(T[j+3]-T[j+1])/2) * (T[j+3]-0) + (((T[j+3]-0)/2.0 - (T[j+2]-T[j])/2.0)/R('wa',(T[j+3]-T[j+1])/2.0 - (T[j+2]-T[j])/2.0)) - (T_int-(T[j+3]-0)/2.0)/R('i',T_int-(T[j+3]-T[j+1])/2) - (T_ext-(T[j+3]-0)/2.0)/R('e',T_ext-(T[j+3]-T[j+1])/2)

   # Derivative of R1 (air balance) with respect to T[j+2] (T[2])
   dR1dT2_a = m_a * Cp('a',(T[j+3]-T[j+1])/2) * (T[j+3]-T[j+1]) + (((T[j+3]-T[j+1])/2.0 - (0-T[j])/2.0)/R('wa',(T[j+3]-T[j+1])/2.0-(T[j+2]-T[j])/2.0)) - (T_int-(T[j+3]-T[j+1])/2.0)/R('i',T_int-(T[j+3]-T[j+1])/2) - (T_ext-(T[j+3]-T[j+1])/2.0)/R('e',T_ext-(T[j+3]-T[j+1])/2)

   # Derivative of R1 (air balance) with respect to T[j+3] (T[3])
   dR1dT3_a = m_a * Cp('a',(T[j+3]-T[j+1])/2) * (0-T[j+1]) + (((0-T[j+1])/2.0 - (T[j+2]-T[j])/2.0)/R('wa',(T[j+3]-T[j+1])/2.0-(T[j+2]-T[j])/2.0)) - (T_int-(0-T[j+1])/2.0)/R('i',T_int-(T[j+3]-T[j+1])/2) - (T_ext-(0-T[j+1])/2.0)/R('e',T_ext-(T[j+3]-T[j+1])/2)

   ##### Region 2: Heat Balance in the Water #####
   # Derivative of R1 (water balance) with respect to T[j] (T[0])
   return dR1dT1_a