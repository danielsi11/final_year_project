function [rate,p_c1,p_c2,p_c3,p_1,p_2,p_3,p_4,p_5,p_6,p_7,p_8,p_9] = ...
RS_perfect_optimisation(H_eff,PAC,g_k,g_ck,W_k,W_ck,P_n_eff,G)
    
[N,K,L] = size(H_eff); 

D1 = zeros(N,N);
D2 = zeros(N,N);
D3 = zeros(N,N);

D1(1,1) = 1;
D2(2,2) = 1;
D3(3,3) = 1;

%CVX Optimisation Tool
cvx_begin quiet

variable p_1(N,1) complex
variable p_2(N,1) complex
variable p_3(N,1) complex
variable p_4(N,1) complex
variable p_5(N,1) complex
variable p_6(N,1) complex
variable p_7(N,1) complex
variable p_8(N,1) complex
variable p_9(N,1) complex
variable p_c1(N,1) complex
variable p_c2(N,1) complex
variable p_c3(N,1) complex
variable C_1 
variable C_2 
variable C_3
variable C_4 
variable C_5 
variable C_6
variable C_7
variable C_8
variable C_9
variable r_1 
variable r_2 
variable r_3
variable r_4 
variable r_5 
variable r_6
variable r_7
variable r_8
variable r_9
variable r_g

expression constraints(1,63);

%Received Power
%T_k Calculation
T_1 = square_abs(H_eff(:,1,2)'*p_c2) + square_abs(H_eff(:,1,3)'*p_c3) + ...
square_abs(H_eff(:,1,1)'*p_1) + square_abs(H_eff(:,1,1)'*p_2) + square_abs(H_eff(:,1,1)'*p_3) + ...
square_abs(H_eff(:,1,2)'*p_4) + square_abs(H_eff(:,1,2)'*p_5) + square_abs(H_eff(:,1,2)'*p_6) + ...
square_abs(H_eff(:,1,3)'*p_7) + square_abs(H_eff(:,1,3)'*p_8) + square_abs(H_eff(:,1,3)'*p_9) + P_n_eff(1);

T_2 = square_abs(H_eff(:,2,2)'*p_c2) + square_abs(H_eff(:,2,3)'*p_c3) + ...
square_abs(H_eff(:,2,1)'*p_1) + square_abs(H_eff(:,2,1)'*p_2) + square_abs(H_eff(:,2,1)'*p_3) + ...
square_abs(H_eff(:,2,2)'*p_4) + square_abs(H_eff(:,2,2)'*p_5) + square_abs(H_eff(:,2,2)'*p_6) + ...
square_abs(H_eff(:,2,3)'*p_7) + square_abs(H_eff(:,2,3)'*p_8) + square_abs(H_eff(:,2,3)'*p_9) + P_n_eff(2);

T_3 = square_abs(H_eff(:,3,2)'*p_c2) + square_abs(H_eff(:,3,3)'*p_c3) + ...
square_abs(H_eff(:,3,1)'*p_2) + square_abs(H_eff(:,3,1)'*p_1) + square_abs(H_eff(:,3,1)'*p_3) + ...
square_abs(H_eff(:,3,2)'*p_4) + square_abs(H_eff(:,3,2)'*p_5) + square_abs(H_eff(:,3,2)'*p_6) + ...
square_abs(H_eff(:,3,3)'*p_7) + square_abs(H_eff(:,3,3)'*p_8) + square_abs(H_eff(:,3,3)'*p_9) + P_n_eff(3);

T_4 = square_abs(H_eff(:,4,2)'*p_c2) + square_abs(H_eff(:,4,3)'*p_c3) + ...
square_abs(H_eff(:,4,1)'*p_2) + square_abs(H_eff(:,4,1)'*p_1) + square_abs(H_eff(:,4,1)'*p_3) + ...
square_abs(H_eff(:,4,2)'*p_4) + square_abs(H_eff(:,4,2)'*p_5) + square_abs(H_eff(:,4,2)'*p_6) + ...
square_abs(H_eff(:,4,3)'*p_7) + square_abs(H_eff(:,4,3)'*p_8) + square_abs(H_eff(:,4,3)'*p_9) + P_n_eff(4);

T_5 = square_abs(H_eff(:,5,2)'*p_c2) + square_abs(H_eff(:,5,3)'*p_c3) + ...
square_abs(H_eff(:,5,1)'*p_3) + square_abs(H_eff(:,5,1)'*p_1) + square_abs(H_eff(:,5,1)'*p_2) + ...
square_abs(H_eff(:,5,2)'*p_4) + square_abs(H_eff(:,5,2)'*p_5) + square_abs(H_eff(:,5,2)'*p_6) + ...
square_abs(H_eff(:,5,3)'*p_7) + square_abs(H_eff(:,5,3)'*p_8) + square_abs(H_eff(:,5,3)'*p_9) + P_n_eff(5);

T_6 = square_abs(H_eff(:,6,2)'*p_c2) + square_abs(H_eff(:,6,3)'*p_c3) + ...
square_abs(H_eff(:,6,1)'*p_3) + square_abs(H_eff(:,6,1)'*p_1) + square_abs(H_eff(:,6,1)'*p_2) + ...
square_abs(H_eff(:,6,2)'*p_4) + square_abs(H_eff(:,6,2)'*p_5) + square_abs(H_eff(:,6,2)'*p_6) + ...
square_abs(H_eff(:,6,3)'*p_7) + square_abs(H_eff(:,6,3)'*p_8) + square_abs(H_eff(:,6,3)'*p_9) + P_n_eff(6);

T_7 = square_abs(H_eff(:,7,1)'*p_c1) + square_abs(H_eff(:,7,3)'*p_c3) + ...
square_abs(H_eff(:,7,2)'*p_4) + square_abs(H_eff(:,7,1)'*p_1) + square_abs(H_eff(:,7,1)'*p_2) + ...
square_abs(H_eff(:,7,1)'*p_3) + square_abs(H_eff(:,7,2)'*p_5) + square_abs(H_eff(:,7,2)'*p_6) + ...
square_abs(H_eff(:,7,3)'*p_7) + square_abs(H_eff(:,7,3)'*p_8) + square_abs(H_eff(:,7,3)'*p_9) + P_n_eff(7);

T_8 = square_abs(H_eff(:,8,1)'*p_c1) + square_abs(H_eff(:,8,3)'*p_c3) + ...
square_abs(H_eff(:,8,2)'*p_4) + square_abs(H_eff(:,8,1)'*p_1) + square_abs(H_eff(:,8,1)'*p_2) + ...
square_abs(H_eff(:,8,1)'*p_3) + square_abs(H_eff(:,8,2)'*p_5) + square_abs(H_eff(:,8,2)'*p_6) + ...
square_abs(H_eff(:,8,3)'*p_7) + square_abs(H_eff(:,8,3)'*p_8) + square_abs(H_eff(:,8,3)'*p_9) + P_n_eff(8);

T_9 = square_abs(H_eff(:,9,1)'*p_c1) + square_abs(H_eff(:,9,3)'*p_c3) + ...
square_abs(H_eff(:,9,2)'*p_5) + square_abs(H_eff(:,9,1)'*p_1) + square_abs(H_eff(:,9,1)'*p_2) + ...
square_abs(H_eff(:,9,1)'*p_3) + square_abs(H_eff(:,9,2)'*p_4) + square_abs(H_eff(:,9,2)'*p_6) + ...
square_abs(H_eff(:,9,3)'*p_7) + square_abs(H_eff(:,9,3)'*p_8) + square_abs(H_eff(:,9,3)'*p_9) + P_n_eff(9);

T_10 = square_abs(H_eff(:,10,1)'*p_c1) + square_abs(H_eff(:,10,3)'*p_c3) + ...
square_abs(H_eff(:,10,2)'*p_5) + square_abs(H_eff(:,10,1)'*p_1) + square_abs(H_eff(:,10,1)'*p_2) + ...
square_abs(H_eff(:,10,1)'*p_3) + square_abs(H_eff(:,10,2)'*p_4) + square_abs(H_eff(:,10,2)'*p_6) + ...
square_abs(H_eff(:,10,3)'*p_7) + square_abs(H_eff(:,10,3)'*p_8) + square_abs(H_eff(:,10,3)'*p_9) + P_n_eff(10);

T_11 = square_abs(H_eff(:,11,1)'*p_c1) + square_abs(H_eff(:,11,3)'*p_c3) + ...
square_abs(H_eff(:,11,2)'*p_6) + square_abs(H_eff(:,11,1)'*p_1) + square_abs(H_eff(:,11,1)'*p_2) + ...
square_abs(H_eff(:,11,1)'*p_3) + square_abs(H_eff(:,11,2)'*p_4) + square_abs(H_eff(:,11,2)'*p_5) + ...
square_abs(H_eff(:,11,3)'*p_7) + square_abs(H_eff(:,11,3)'*p_8) + square_abs(H_eff(:,11,3)'*p_9) + P_n_eff(11);

T_12 = square_abs(H_eff(:,12,1)'*p_c1) + square_abs(H_eff(:,12,3)'*p_c3) + ...
square_abs(H_eff(:,12,2)'*p_6) + square_abs(H_eff(:,12,1)'*p_1) + square_abs(H_eff(:,12,1)'*p_2) + ...
square_abs(H_eff(:,12,1)'*p_3) + square_abs(H_eff(:,12,2)'*p_4) + square_abs(H_eff(:,12,2)'*p_5) + ...
square_abs(H_eff(:,12,3)'*p_7) + square_abs(H_eff(:,12,3)'*p_8) + square_abs(H_eff(:,12,3)'*p_9) + P_n_eff(12);

T_13 = square_abs(H_eff(:,13,1)'*p_c1) + square_abs(H_eff(:,13,2)'*p_c2) + ...
square_abs(H_eff(:,13,3)'*p_7) + square_abs(H_eff(:,13,1)'*p_1) + square_abs(H_eff(:,13,1)'*p_2) + ...
square_abs(H_eff(:,13,1)'*p_3) + square_abs(H_eff(:,13,2)'*p_4) + square_abs(H_eff(:,13,2)'*p_5) + ...
square_abs(H_eff(:,13,2)'*p_6) + square_abs(H_eff(:,13,3)'*p_8) + square_abs(H_eff(:,13,3)'*p_9) + P_n_eff(13);

T_14 = square_abs(H_eff(:,14,1)'*p_c1) + square_abs(H_eff(:,14,2)'*p_c2) + ...
square_abs(H_eff(:,14,3)'*p_7) + square_abs(H_eff(:,14,1)'*p_1) + square_abs(H_eff(:,14,1)'*p_2) + ...
square_abs(H_eff(:,14,1)'*p_3) + square_abs(H_eff(:,14,2)'*p_4) + square_abs(H_eff(:,14,2)'*p_5) + ...
square_abs(H_eff(:,14,2)'*p_6) + square_abs(H_eff(:,14,3)'*p_8) + square_abs(H_eff(:,14,3)'*p_9) + P_n_eff(14);

T_15 = square_abs(H_eff(:,15,1)'*p_c1) + square_abs(H_eff(:,15,2)'*p_c2) + ...
square_abs(H_eff(:,15,3)'*p_8) + square_abs(H_eff(:,15,1)'*p_1) + square_abs(H_eff(:,15,1)'*p_2) + ...
square_abs(H_eff(:,15,1)'*p_3) + square_abs(H_eff(:,15,2)'*p_4) + square_abs(H_eff(:,15,2)'*p_5) + ...
square_abs(H_eff(:,15,2)'*p_6) + square_abs(H_eff(:,15,3)'*p_7) + square_abs(H_eff(:,15,3)'*p_9) + P_n_eff(15);

T_16 = square_abs(H_eff(:,16,1)'*p_c1) + square_abs(H_eff(:,16,2)'*p_c2) + ...
square_abs(H_eff(:,16,3)'*p_8) + square_abs(H_eff(:,16,1)'*p_1) + square_abs(H_eff(:,16,1)'*p_2) + ...
square_abs(H_eff(:,16,1)'*p_3) + square_abs(H_eff(:,16,2)'*p_4) + square_abs(H_eff(:,16,2)'*p_5) + ...
square_abs(H_eff(:,16,2)'*p_6) + square_abs(H_eff(:,16,3)'*p_7) + square_abs(H_eff(:,16,3)'*p_9) + P_n_eff(16);

T_17 = square_abs(H_eff(:,17,1)'*p_c1) + square_abs(H_eff(:,17,2)'*p_c2) + ...
square_abs(H_eff(:,17,3)'*p_9) + square_abs(H_eff(:,17,1)'*p_1) + square_abs(H_eff(:,17,1)'*p_2) + ...
square_abs(H_eff(:,17,1)'*p_3) + square_abs(H_eff(:,17,2)'*p_4) + square_abs(H_eff(:,17,2)'*p_5) + ...
square_abs(H_eff(:,17,2)'*p_6) + square_abs(H_eff(:,17,3)'*p_7) + square_abs(H_eff(:,17,3)'*p_8) + P_n_eff(17);

T_18 = square_abs(H_eff(:,18,1)'*p_c1) + square_abs(H_eff(:,18,2)'*p_c2) + ...
square_abs(H_eff(:,18,3)'*p_9) + square_abs(H_eff(:,18,1)'*p_1) + square_abs(H_eff(:,18,1)'*p_2) + ...
square_abs(H_eff(:,18,1)'*p_3) + square_abs(H_eff(:,18,2)'*p_4) + square_abs(H_eff(:,18,2)'*p_5) + ...
square_abs(H_eff(:,18,2)'*p_6) + square_abs(H_eff(:,18,3)'*p_7) + square_abs(H_eff(:,18,3)'*p_8) + P_n_eff(18);

%T_ck Calculation
T_c1 = square_abs(H_eff(:,1,1)'*p_c1) + T_1;
T_c2 = square_abs(H_eff(:,2,1)'*p_c1) + T_2;
T_c3 = square_abs(H_eff(:,3,1)'*p_c1) + T_3;
T_c4 = square_abs(H_eff(:,4,1)'*p_c1) + T_4;
T_c5 = square_abs(H_eff(:,5,1)'*p_c1) + T_5;
T_c6 = square_abs(H_eff(:,6,1)'*p_c1) + T_6;
T_c7 = square_abs(H_eff(:,7,2)'*p_c2) + T_7;
T_c8 = square_abs(H_eff(:,8,2)'*p_c2) + T_8;
T_c9 = square_abs(H_eff(:,9,2)'*p_c2) + T_9;
T_c10 = square_abs(H_eff(:,10,2)'*p_c2) + T_10;
T_c11 = square_abs(H_eff(:,11,2)'*p_c2) + T_11;
T_c12 = square_abs(H_eff(:,12,2)'*p_c2) + T_12;
T_c13 = square_abs(H_eff(:,13,3)'*p_c3) + T_13;
T_c14 = square_abs(H_eff(:,14,3)'*p_c3) + T_14;
T_c15 = square_abs(H_eff(:,15,3)'*p_c3) + T_15;
T_c16 = square_abs(H_eff(:,16,3)'*p_c3) + T_16;
T_c17 = square_abs(H_eff(:,17,3)'*p_c3) + T_17;
T_c18 = square_abs(H_eff(:,18,3)'*p_c3) + T_18;

%MSE
%E_k Calculation
E_1 = abs(g_k(1))^2*T_1 - 2*real(g_k(1)*H_eff(:,1,1)'*p_1) + 1;
E_2 = abs(g_k(2))^2*T_2 - 2*real(g_k(2)*H_eff(:,2,1)'*p_1) + 1;
E_3 = abs(g_k(3))^2*T_3 - 2*real(g_k(3)*H_eff(:,3,1)'*p_2) + 1;
E_4 = abs(g_k(4))^2*T_4 - 2*real(g_k(4)*H_eff(:,4,1)'*p_2) + 1;
E_5 = abs(g_k(5))^2*T_5 - 2*real(g_k(5)*H_eff(:,5,1)'*p_3) + 1;
E_6 = abs(g_k(6))^2*T_6 - 2*real(g_k(6)*H_eff(:,6,1)'*p_3) + 1;
E_7 = abs(g_k(7))^2*T_7 - 2*real(g_k(7)*H_eff(:,7,2)'*p_4) + 1;
E_8 = abs(g_k(8))^2*T_8 - 2*real(g_k(8)*H_eff(:,8,2)'*p_4) + 1;
E_9 = abs(g_k(9))^2*T_9 - 2*real(g_k(9)*H_eff(:,9,2)'*p_5) + 1;
E_10 = abs(g_k(10))^2*T_10 - 2*real(g_k(10)*H_eff(:,10,2)'*p_5) + 1;
E_11 = abs(g_k(11))^2*T_11 - 2*real(g_k(11)*H_eff(:,11,2)'*p_6) + 1;
E_12 = abs(g_k(12))^2*T_12 - 2*real(g_k(12)*H_eff(:,12,2)'*p_6) + 1;
E_13 = abs(g_k(13))^2*T_13 - 2*real(g_k(13)*H_eff(:,13,3)'*p_7) + 1;
E_14 = abs(g_k(14))^2*T_14 - 2*real(g_k(14)*H_eff(:,14,3)'*p_7) + 1;
E_15 = abs(g_k(15))^2*T_15 - 2*real(g_k(15)*H_eff(:,15,3)'*p_8) + 1;
E_16 = abs(g_k(16))^2*T_16 - 2*real(g_k(16)*H_eff(:,16,3)'*p_8) + 1;
E_17 = abs(g_k(17))^2*T_17 - 2*real(g_k(17)*H_eff(:,17,3)'*p_9) + 1;
E_18 = abs(g_k(18))^2*T_18 - 2*real(g_k(18)*H_eff(:,18,3)'*p_9) + 1;

%E_ck Calculation
E_c1 = abs(g_ck(1))^2*T_c1 - 2*real(g_ck(1)*H_eff(:,1,1)'*p_c1) + 1;
E_c2 = abs(g_ck(2))^2*T_c2 - 2*real(g_ck(2)*H_eff(:,2,1)'*p_c1) + 1;
E_c3 = abs(g_ck(3))^2*T_c3 - 2*real(g_ck(3)*H_eff(:,3,1)'*p_c1) + 1;
E_c4 = abs(g_ck(4))^2*T_c4 - 2*real(g_ck(4)*H_eff(:,4,1)'*p_c1) + 1;
E_c5 = abs(g_ck(5))^2*T_c5 - 2*real(g_ck(5)*H_eff(:,5,1)'*p_c1) + 1;
E_c6 = abs(g_ck(6))^2*T_c6 - 2*real(g_ck(6)*H_eff(:,6,1)'*p_c1) + 1;
E_c7 = abs(g_ck(7))^2*T_c7 - 2*real(g_ck(7)*H_eff(:,7,2)'*p_c2) + 1;
E_c8 = abs(g_ck(8))^2*T_c8 - 2*real(g_ck(8)*H_eff(:,8,2)'*p_c2) + 1;
E_c9 = abs(g_ck(9))^2*T_c9 - 2*real(g_ck(9)*H_eff(:,9,2)'*p_c2) + 1;
E_c10 = abs(g_ck(10))^2*T_c10 - 2*real(g_ck(10)*H_eff(:,10,2)'*p_c2) + 1;
E_c11 = abs(g_ck(11))^2*T_c11 - 2*real(g_ck(11)*H_eff(:,11,2)'*p_c2) + 1;
E_c12 = abs(g_ck(12))^2*T_c12 - 2*real(g_ck(12)*H_eff(:,12,2)'*p_c2) + 1;
E_c13 = abs(g_ck(13))^2*T_c13 - 2*real(g_ck(13)*H_eff(:,13,3)'*p_c3) + 1;
E_c14 = abs(g_ck(14))^2*T_c14 - 2*real(g_ck(14)*H_eff(:,14,3)'*p_c3) + 1;
E_c15 = abs(g_ck(15))^2*T_c15 - 2*real(g_ck(15)*H_eff(:,15,3)'*p_c3) + 1;
E_c16 = abs(g_ck(16))^2*T_c16 - 2*real(g_ck(16)*H_eff(:,16,3)'*p_c3) + 1;
E_c17 = abs(g_ck(17))^2*T_c17 - 2*real(g_ck(17)*H_eff(:,17,3)'*p_c3) + 1;
E_c18 = abs(g_ck(18))^2*T_c18 - 2*real(g_ck(18)*H_eff(:,18,3)'*p_c3) + 1;

%Rate-WMMSE Relationship
%Private Stream
X_1 = W_k(1)*E_1 - log2(W_k(1));
X_2 = W_k(2)*E_2 - log2(W_k(2));
X_3 = W_k(3)*E_3 - log2(W_k(3));
X_4 = W_k(4)*E_4 - log2(W_k(4));
X_5 = W_k(5)*E_5 - log2(W_k(5));
X_6 = W_k(6)*E_6 - log2(W_k(6));
X_7 = W_k(7)*E_7 - log2(W_k(7));
X_8 = W_k(8)*E_8 - log2(W_k(8));
X_9 = W_k(9)*E_9 - log2(W_k(9));
X_10 = W_k(10)*E_10 - log2(W_k(10));
X_11 = W_k(11)*E_11 - log2(W_k(11));
X_12 = W_k(12)*E_12 - log2(W_k(12));
X_13 = W_k(13)*E_13 - log2(W_k(13));
X_14 = W_k(14)*E_14 - log2(W_k(14));
X_15 = W_k(15)*E_15 - log2(W_k(15));
X_16 = W_k(16)*E_16 - log2(W_k(16));
X_17 = W_k(17)*E_17 - log2(W_k(17));
X_18 = W_k(18)*E_18 - log2(W_k(18));

%Common Stream
X_c1 = W_ck(1)*E_c1 - log2(W_ck(1));
X_c2 = W_ck(2)*E_c2 - log2(W_ck(2));
X_c3 = W_ck(3)*E_c3 - log2(W_ck(3));
X_c4 = W_ck(4)*E_c4 - log2(W_ck(4));
X_c5 = W_ck(5)*E_c5 - log2(W_ck(5));
X_c6 = W_ck(6)*E_c6 - log2(W_ck(6));
X_c7 = W_ck(7)*E_c7 - log2(W_ck(7));
X_c8 = W_ck(8)*E_c8 - log2(W_ck(8));
X_c9 = W_ck(9)*E_c9 - log2(W_ck(9));
X_c10 = W_ck(10)*E_c10 - log2(W_ck(10));
X_c11 = W_ck(11)*E_c11 - log2(W_ck(11));
X_c12 = W_ck(12)*E_c12 - log2(W_ck(12));
X_c13 = W_ck(13)*E_c13 - log2(W_ck(13));
X_c14 = W_ck(14)*E_c14 - log2(W_ck(14));
X_c15 = W_ck(15)*E_c15 - log2(W_ck(15));
X_c16 = W_ck(16)*E_c16 - log2(W_ck(16));
X_c17 = W_ck(17)*E_c17 - log2(W_ck(17));
X_c18 = W_ck(18)*E_c18 - log2(W_ck(18));

%Objective Function
object_func = r_g;

%Optimisation
maximize(object_func)

%Constraints    
constraints(1) = C_1 + r_1 - r_g;
constraints(2) = C_2 + r_2 - r_g;
constraints(3) = C_3 + r_3 - r_g;
constraints(4) = C_4 + r_4 - r_g;
constraints(5) = C_5 + r_5 - r_g;
constraints(6) = C_6 + r_6 - r_g;
constraints(7) = C_7 + r_7 - r_g;
constraints(8) = C_8 + r_8 - r_g;
constraints(9) = C_9 + r_9 - r_g;
constraints(10) = 1 - X_1 - r_1;
constraints(11) = 1 - X_2 - r_1;
constraints(12) = 1 - X_3 - r_2;
constraints(13) = 1 - X_4 - r_2;
constraints(14) = 1 - X_5 - r_3;
constraints(15) = 1 - X_6 - r_3;
constraints(16) = 1 - X_7 - r_4;
constraints(17) = 1 - X_8 - r_4;
constraints(18) = 1 - X_9 - r_5;
constraints(19) = 1 - X_10 - r_5;
constraints(20) = 1 - X_11 - r_6;
constraints(21) = 1 - X_12 - r_6;
constraints(22) = 1 - X_13 - r_7;
constraints(23) = 1 - X_14 - r_7;
constraints(24) = 1 - X_15 - r_8;
constraints(25) = 1 - X_16 - r_8;
constraints(26) = 1 - X_17 - r_9;
constraints(27) = 1 - X_18 - r_9;
constraints(28) = 1 - X_c1 - C_1 - C_2 - C_3;
constraints(29) = 1 - X_c2 - C_1 - C_2 - C_3;
constraints(30) = 1 - X_c3 - C_1 - C_2 - C_3;
constraints(31) = 1 - X_c4 - C_1 - C_2 - C_3;
constraints(32) = 1 - X_c5 - C_1 - C_2 - C_3;
constraints(33) = 1 - X_c6 - C_1 - C_2 - C_3;
constraints(34) = 1 - X_c7 - C_4 - C_5 - C_6;
constraints(35) = 1 - X_c8 - C_4 - C_5 - C_6;
constraints(36) = 1 - X_c9 - C_4 - C_5 - C_6;
constraints(37) = 1 - X_c10 - C_4 - C_5 - C_6;
constraints(38) = 1 - X_c11 - C_4 - C_5 - C_6;
constraints(39) = 1 - X_c12 - C_4 - C_5 - C_6;
constraints(40) = 1 - X_c13 - C_7 - C_8 - C_9;
constraints(41) = 1 - X_c14 - C_7 - C_8 - C_9;
constraints(42) = 1 - X_c15 - C_7 - C_8 - C_9;
constraints(43) = 1 - X_c16 - C_7 - C_8 - C_9;
constraints(44) = 1 - X_c17 - C_7 - C_8 - C_9;
constraints(45) = 1 - X_c18 - C_7 - C_8 - C_9;
constraints(46) = C_1;
constraints(47) = C_2;
constraints(48) = C_3;
constraints(49) = C_4;
constraints(50) = C_5;
constraints(51) = C_6;
constraints(52) = C_7;
constraints(53) = C_8;
constraints(54) = C_9;

constraints(55) = PAC - (p_c1'*G(:,:,1,1)'*D1*G(:,:,1,1)*p_c1 + p_1'*G(:,:,1,1)'*D1*G(:,:,1,1)*p_1 + p_2'*G(:,:,1,1)'*D1*G(:,:,1,1)*p_2 + p_3'*G(:,:,1,1)'*D1*G(:,:,1,1)*p_3 + ... 
p_c2'*G(:,:,1,2)'*D1*G(:,:,1,2)*p_c2 + p_4'*G(:,:,1,2)'*D1*G(:,:,1,2)*p_4 + p_5'*G(:,:,1,2)'*D1*G(:,:,1,2)*p_5 + p_6'*G(:,:,1,2)'*D1*G(:,:,1,2)*p_6 + ...
p_c3'*G(:,:,1,3)'*D1*G(:,:,1,3)*p_c3 + p_7'*G(:,:,1,3)'*D1*G(:,:,1,3)*p_7 + p_8'*G(:,:,1,3)'*D1*G(:,:,1,3)*p_8 + p_9'*G(:,:,1,3)'*D1*G(:,:,1,3)*p_9 + 1);

constraints(56) = PAC - (p_c1'*G(:,:,1,1)'*D2*G(:,:,1,1)*p_c1 + p_1'*G(:,:,1,1)'*D2*G(:,:,1,1)*p_1 + p_2'*G(:,:,1,1)'*D2*G(:,:,1,1)*p_2 + p_3'*G(:,:,1,1)'*D2*G(:,:,1,1)*p_3 + ... 
p_c2'*G(:,:,1,2)'*D2*G(:,:,1,2)*p_c2 + p_4'*G(:,:,1,2)'*D2*G(:,:,1,2)*p_4 + p_5'*G(:,:,1,2)'*D2*G(:,:,1,2)*p_5 + p_6'*G(:,:,1,2)'*D2*G(:,:,1,2)*p_6 + ...
p_c3'*G(:,:,1,3)'*D2*G(:,:,1,3)*p_c3 + p_7'*G(:,:,1,3)'*D2*G(:,:,1,3)*p_7 + p_8'*G(:,:,1,3)'*D2*G(:,:,1,3)*p_8 + p_9'*G(:,:,1,3)'*D2*G(:,:,1,3)*p_9 + 1);

constraints(57) = PAC - (p_c1'*G(:,:,1,1)'*D3*G(:,:,1,1)*p_c1 + p_1'*G(:,:,1,1)'*D3*G(:,:,1,1)*p_1 + p_2'*G(:,:,1,1)'*D3*G(:,:,1,1)*p_2 + p_3'*G(:,:,1,1)'*D3*G(:,:,1,1)*p_3 + ... 
p_c2'*G(:,:,1,2)'*D3*G(:,:,1,2)*p_c2 + p_4'*G(:,:,1,2)'*D3*G(:,:,1,2)*p_4 + p_5'*G(:,:,1,2)'*D3*G(:,:,1,2)*p_5 + p_6'*G(:,:,1,2)'*D3*G(:,:,1,2)*p_6 + ...
p_c3'*G(:,:,1,3)'*D3*G(:,:,1,3)*p_c3 + p_7'*G(:,:,1,3)'*D3*G(:,:,1,3)*p_7 + p_8'*G(:,:,1,3)'*D3*G(:,:,1,3)*p_8 + p_9'*G(:,:,1,3)'*D3*G(:,:,1,3)*p_9 + 1);

constraints(58) = PAC - (p_c1'*G(:,:,2,1)'*D1*G(:,:,2,1)*p_c1 + p_1'*G(:,:,2,1)'*D1*G(:,:,2,1)*p_1 + p_2'*G(:,:,2,1)'*D1*G(:,:,2,1)*p_2 + p_3'*G(:,:,2,1)'*D1*G(:,:,2,1)*p_3 + ... 
p_c2'*G(:,:,2,2)'*D1*G(:,:,2,2)*p_c2 + p_4'*G(:,:,2,2)'*D1*G(:,:,2,2)*p_4 + p_5'*G(:,:,2,2)'*D1*G(:,:,2,2)*p_5 + p_6'*G(:,:,2,2)'*D1*G(:,:,2,2)*p_6 + ...
p_c3'*G(:,:,2,3)'*D1*G(:,:,2,3)*p_c3 + p_7'*G(:,:,2,3)'*D1*G(:,:,2,3)*p_7 + p_8'*G(:,:,2,3)'*D1*G(:,:,2,3)*p_8 + p_9'*G(:,:,2,3)'*D1*G(:,:,2,3)*p_9 + 1);

constraints(59) = PAC - (p_c1'*G(:,:,2,1)'*D2*G(:,:,2,1)*p_c1 + p_1'*G(:,:,2,1)'*D2*G(:,:,2,1)*p_1 + p_2'*G(:,:,2,1)'*D2*G(:,:,2,1)*p_2 + p_3'*G(:,:,2,1)'*D2*G(:,:,2,1)*p_3 + ... 
p_c2'*G(:,:,2,2)'*D2*G(:,:,2,2)*p_c2 + p_4'*G(:,:,2,2)'*D2*G(:,:,2,2)*p_4 + p_5'*G(:,:,2,2)'*D2*G(:,:,2,2)*p_5 + p_6'*G(:,:,2,2)'*D2*G(:,:,2,2)*p_6 + ...
p_c3'*G(:,:,2,3)'*D2*G(:,:,2,3)*p_c3 + p_7'*G(:,:,2,3)'*D2*G(:,:,2,3)*p_7 + p_8'*G(:,:,2,3)'*D2*G(:,:,2,3)*p_8 + p_9'*G(:,:,2,3)'*D2*G(:,:,2,3)*p_9 + 1);

constraints(60) = PAC - (p_c1'*G(:,:,2,1)'*D3*G(:,:,2,1)*p_c1 + p_1'*G(:,:,2,1)'*D3*G(:,:,2,1)*p_1 + p_2'*G(:,:,2,1)'*D3*G(:,:,2,1)*p_2 + p_3'*G(:,:,2,1)'*D3*G(:,:,2,1)*p_3 + ... 
p_c2'*G(:,:,2,2)'*D3*G(:,:,2,2)*p_c2 + p_4'*G(:,:,2,2)'*D3*G(:,:,2,2)*p_4 + p_5'*G(:,:,2,2)'*D3*G(:,:,2,2)*p_5 + p_6'*G(:,:,2,2)'*D3*G(:,:,2,2)*p_6 + ...
p_c3'*G(:,:,2,3)'*D3*G(:,:,2,3)*p_c3 + p_7'*G(:,:,2,3)'*D3*G(:,:,2,3)*p_7 + p_8'*G(:,:,2,3)'*D3*G(:,:,2,3)*p_8 + p_9'*G(:,:,2,3)'*D3*G(:,:,2,3)*p_9 + 1);

constraints(61) = PAC - (p_c1'*G(:,:,3,1)'*D1*G(:,:,3,1)*p_c1 + p_1'*G(:,:,3,1)'*D1*G(:,:,3,1)*p_1 + p_2'*G(:,:,3,1)'*D1*G(:,:,3,1)*p_2 + p_3'*G(:,:,3,1)'*D1*G(:,:,3,1)*p_3 + ... 
p_c2'*G(:,:,3,2)'*D1*G(:,:,3,2)*p_c2 + p_4'*G(:,:,3,2)'*D1*G(:,:,3,2)*p_4 + p_5'*G(:,:,3,2)'*D1*G(:,:,3,2)*p_5 + p_6'*G(:,:,3,2)'*D1*G(:,:,3,2)*p_6 + ...
p_c3'*G(:,:,3,3)'*D1*G(:,:,3,3)*p_c3 + p_7'*G(:,:,3,3)'*D1*G(:,:,3,3)*p_7 + p_8'*G(:,:,3,3)'*D1*G(:,:,3,3)*p_8 + p_9'*G(:,:,3,3)'*D1*G(:,:,3,3)*p_9 + 1);

constraints(62) = PAC - (p_c1'*G(:,:,3,1)'*D2*G(:,:,3,1)*p_c1 + p_1'*G(:,:,3,1)'*D2*G(:,:,3,1)*p_1 + p_2'*G(:,:,3,1)'*D2*G(:,:,3,1)*p_2 + p_3'*G(:,:,3,1)'*D2*G(:,:,3,1)*p_3 + ... 
p_c2'*G(:,:,3,2)'*D2*G(:,:,3,2)*p_c2 + p_4'*G(:,:,3,2)'*D2*G(:,:,3,2)*p_4 + p_5'*G(:,:,3,2)'*D2*G(:,:,3,2)*p_5 + p_6'*G(:,:,3,2)'*D2*G(:,:,3,2)*p_6 + ...
p_c3'*G(:,:,3,3)'*D2*G(:,:,3,3)*p_c3 + p_7'*G(:,:,3,3)'*D2*G(:,:,3,3)*p_7 + p_8'*G(:,:,3,3)'*D2*G(:,:,3,3)*p_8 + p_9'*G(:,:,3,3)'*D2*G(:,:,3,3)*p_9 + 1);

constraints(63) = PAC - (p_c1'*G(:,:,3,1)'*D3*G(:,:,3,1)*p_c1 + p_1'*G(:,:,3,1)'*D3*G(:,:,3,1)*p_1 + p_2'*G(:,:,3,1)'*D3*G(:,:,3,1)*p_2 + p_3'*G(:,:,3,1)'*D3*G(:,:,3,1)*p_3 + ... 
p_c2'*G(:,:,3,2)'*D3*G(:,:,3,2)*p_c2 + p_4'*G(:,:,3,2)'*D3*G(:,:,3,2)*p_4 + p_5'*G(:,:,3,2)'*D3*G(:,:,3,2)*p_5 + p_6'*G(:,:,3,2)'*D3*G(:,:,3,2)*p_6 + ...
p_c3'*G(:,:,3,3)'*D3*G(:,:,3,3)*p_c3 + p_7'*G(:,:,3,3)'*D3*G(:,:,3,3)*p_7 + p_8'*G(:,:,3,3)'*D3*G(:,:,3,3)*p_8 + p_9'*G(:,:,3,3)'*D3*G(:,:,3,3)*p_9 + 1);


subject to
    constraints >= zeros(size(constraints))

cvx_end

%Result
rate = object_func;
