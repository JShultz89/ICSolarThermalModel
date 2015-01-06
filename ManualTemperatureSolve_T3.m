clear;
clc;
format long;

T_0 = 13;
T_1 = 20;
T_int = 22.5;
T_ext = 25;
m_w = 0.00084931862198712224;
m_a = 0.3597862499999999;
Cp_w = 4.188774760737728;
Cp_a = 1.005; 
R_pipe = 1472.0223510771341;
R_int = 0.52972312781694775;
R_ext = 0.10670725480107474;

A_12 = m_w * Cp_w + 1/2/R_pipe;

A_13 = -1/2/R_pipe; 

F_1 = T_0 * (m_w * Cp_w - 1/2/R_pipe) ...
    + T_1 * 1/2/R_pipe;


A_22 = -1/2/R_pipe;

A_23 = m_a * Cp_a + 1/2/R_pipe + 1/R_int + 1/R_ext;

F_2 = T_1 * ( m_a * Cp_a - 1.0/2/R_pipe) ...
    + T_0 * (1/2/R_pipe) ...
    + T_int * (1/R_int) ...
    + T_ext * (1/R_ext);

% A * T = F
% T = [T_2; T_3]

A = [A_12 A_13; A_22 A_23]

F = [F_1; F_2]

T = A \ F;

T_2 = T(1)
T_3 = T(2)
