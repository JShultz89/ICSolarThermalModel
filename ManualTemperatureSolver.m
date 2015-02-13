% clear;
% clc;
format long;

T_0 = 13;
T_1 = 20;
T_int = 22.5;
T_ext = 25;
m_w = 8.5*10^(-7) * 998.8376 % mass flow rate = volumetric flowrate * density
m_a = 2.0 * 1.1992875 * 0.3*0.3 % mass flow rate = velocity * density * cross-section
Cp_w = 4.18759087915;
Cp_a = 1.005; 
% radius_inside = 3*10^(-3);
% length_tubing = 3.0;
% h_w = 719.797570611;
R_pipe = 20.7307482776; % = 1/(2*pi*radius_inside*length_tubing*h_w)
R_int = 5.32864511015;
R_ext = 1.10123830038;

A_12 = m_w * Cp_w + 1/R_pipe;

A_13 = -1/R_pipe; 

F_1 = T_0 * (m_w * Cp_w) ...
    + T_1 * 0;


A_22 = -1/R_pipe;

A_23 = m_a * Cp_a + 1/R_pipe + 1/R_int + 1/R_ext;

F_2 = T_1 * (m_a * Cp_a) ...
    + T_0 * (0) ...
    + T_int * (1/R_int) ...
    + T_ext * (1/R_ext);

% A * T = F
% T = [T_2; T_3]

A = [A_12 A_13; A_22 A_23]

F = [F_1; F_2]

T = A \ F;

T_2 = T(1)
T_3 = T(2)
