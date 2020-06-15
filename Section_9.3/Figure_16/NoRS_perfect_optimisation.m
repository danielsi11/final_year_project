function [rate,v_1,v_2,v_3,v_4,v_5,v_6,v_7,v_8,v_9] = NoRS_perfect_optimisation(H_eff,PAC,g_k,W_k,P_n_eff,G_eff,R)
    
[N,K,L] = size(H_eff); 

D1 = zeros(N,N);
D2 = zeros(N,N);
D3 = zeros(N,N);

D1(1,1) = 1;
D2(2,2) = 1;
D3(3,3) = 1;

%CVX Optimisation Tool
cvx_begin quiet

variable v_1(N,1) complex
variable v_2(N,1) complex
variable v_3(N,1) complex
variable v_4(N,1) complex
variable v_5(N,1) complex
variable v_6(N,1) complex
variable v_7(N,1) complex
variable v_8(N,1) complex
variable v_9(N,1) complex
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

expression constraints(1,36);

%Received Power
%T_k Calculation
T_1 = square_abs(H_eff(:,1,1)'*v_1) + square_abs(H_eff(:,1,1)'*v_2) + square_abs(H_eff(:,1,1)'*v_3) + ...
square_abs(H_eff(:,1,2)'*v_4) + square_abs(H_eff(:,1,2)'*v_5) + square_abs(H_eff(:,1,2)'*v_6) + ...
square_abs(H_eff(:,1,3)'*v_7) + square_abs(H_eff(:,1,3)'*v_8) + square_abs(H_eff(:,1,3)'*v_9) + P_n_eff(1);

T_2 = square_abs(H_eff(:,2,1)'*v_1) + square_abs(H_eff(:,2,1)'*v_2) + square_abs(H_eff(:,2,1)'*v_3) + ...
square_abs(H_eff(:,2,2)'*v_4) + square_abs(H_eff(:,2,2)'*v_5) + square_abs(H_eff(:,2,2)'*v_6) + ...
square_abs(H_eff(:,2,3)'*v_7) + square_abs(H_eff(:,2,3)'*v_8) + square_abs(H_eff(:,2,3)'*v_9) + P_n_eff(2);

T_3 = square_abs(H_eff(:,3,1)'*v_2) + square_abs(H_eff(:,3,1)'*v_1) + square_abs(H_eff(:,3,1)'*v_3) + ...
square_abs(H_eff(:,3,2)'*v_4) + square_abs(H_eff(:,3,2)'*v_5) + square_abs(H_eff(:,3,2)'*v_6) + ...
square_abs(H_eff(:,3,3)'*v_7) + square_abs(H_eff(:,3,3)'*v_8) + square_abs(H_eff(:,3,3)'*v_9) + P_n_eff(3);

T_4 = square_abs(H_eff(:,4,1)'*v_2) + square_abs(H_eff(:,4,1)'*v_1) + square_abs(H_eff(:,4,1)'*v_3) + ...
square_abs(H_eff(:,4,2)'*v_4) + square_abs(H_eff(:,4,2)'*v_5) + square_abs(H_eff(:,4,2)'*v_6) + ...
square_abs(H_eff(:,4,3)'*v_7) + square_abs(H_eff(:,4,3)'*v_8) + square_abs(H_eff(:,4,3)'*v_9) + P_n_eff(4);

T_5 = square_abs(H_eff(:,5,1)'*v_3) + square_abs(H_eff(:,5,1)'*v_1) + square_abs(H_eff(:,5,1)'*v_2) + ...
square_abs(H_eff(:,5,2)'*v_4) + square_abs(H_eff(:,5,2)'*v_5) + square_abs(H_eff(:,5,2)'*v_6) + ...
square_abs(H_eff(:,5,3)'*v_7) + square_abs(H_eff(:,5,3)'*v_8) + square_abs(H_eff(:,5,3)'*v_9) + P_n_eff(5);

T_6 = square_abs(H_eff(:,6,1)'*v_3) + square_abs(H_eff(:,6,1)'*v_1) + square_abs(H_eff(:,6,1)'*v_2) + ...
square_abs(H_eff(:,6,2)'*v_4) + square_abs(H_eff(:,6,2)'*v_5) + square_abs(H_eff(:,6,2)'*v_6) + ...
square_abs(H_eff(:,6,3)'*v_7) + square_abs(H_eff(:,6,3)'*v_8) + square_abs(H_eff(:,6,3)'*v_9) + P_n_eff(6);

T_7 = square_abs(H_eff(:,7,2)'*v_4) + square_abs(H_eff(:,7,1)'*v_1) + square_abs(H_eff(:,7,1)'*v_2) + ...
square_abs(H_eff(:,7,1)'*v_3) + square_abs(H_eff(:,7,2)'*v_5) + square_abs(H_eff(:,7,2)'*v_6) + ...
square_abs(H_eff(:,7,3)'*v_7) + square_abs(H_eff(:,7,3)'*v_8) + square_abs(H_eff(:,7,3)'*v_9) + P_n_eff(7);

T_8 = square_abs(H_eff(:,8,2)'*v_4) + square_abs(H_eff(:,8,1)'*v_1) + square_abs(H_eff(:,8,1)'*v_2) + ...
square_abs(H_eff(:,8,1)'*v_3) + square_abs(H_eff(:,8,2)'*v_5) + square_abs(H_eff(:,8,2)'*v_6) + ...
square_abs(H_eff(:,8,3)'*v_7) + square_abs(H_eff(:,8,3)'*v_8) + square_abs(H_eff(:,8,3)'*v_9) + P_n_eff(8);

T_9 = square_abs(H_eff(:,9,2)'*v_5) + square_abs(H_eff(:,9,1)'*v_1) + square_abs(H_eff(:,9,1)'*v_2) + ...
square_abs(H_eff(:,9,1)'*v_3) + square_abs(H_eff(:,9,2)'*v_4) + square_abs(H_eff(:,9,2)'*v_6) + ...
square_abs(H_eff(:,9,3)'*v_7) + square_abs(H_eff(:,9,3)'*v_8) + square_abs(H_eff(:,9,3)'*v_9) + P_n_eff(9);

T_10 = square_abs(H_eff(:,10,2)'*v_5) + square_abs(H_eff(:,10,1)'*v_1) + square_abs(H_eff(:,10,1)'*v_2) + ...
square_abs(H_eff(:,10,1)'*v_3) + square_abs(H_eff(:,10,2)'*v_4) + square_abs(H_eff(:,10,2)'*v_6) + ...
square_abs(H_eff(:,10,3)'*v_7) + square_abs(H_eff(:,10,3)'*v_8) + square_abs(H_eff(:,10,3)'*v_9) + P_n_eff(10);

T_11 = square_abs(H_eff(:,11,2)'*v_6) + square_abs(H_eff(:,11,1)'*v_1) + square_abs(H_eff(:,11,1)'*v_2) + ...
square_abs(H_eff(:,11,1)'*v_3) + square_abs(H_eff(:,11,2)'*v_4) + square_abs(H_eff(:,11,2)'*v_5) + ...
square_abs(H_eff(:,11,3)'*v_7) + square_abs(H_eff(:,11,3)'*v_8) + square_abs(H_eff(:,11,3)'*v_9) + P_n_eff(11);

T_12 = square_abs(H_eff(:,12,2)'*v_6) + square_abs(H_eff(:,12,1)'*v_1) + square_abs(H_eff(:,12,1)'*v_2) + ...
square_abs(H_eff(:,12,1)'*v_3) + square_abs(H_eff(:,12,2)'*v_4) + square_abs(H_eff(:,12,2)'*v_5) + ...
square_abs(H_eff(:,12,3)'*v_7) + square_abs(H_eff(:,12,3)'*v_8) + square_abs(H_eff(:,12,3)'*v_9) + P_n_eff(12);

T_13 = square_abs(H_eff(:,13,3)'*v_7) + square_abs(H_eff(:,13,1)'*v_1) + square_abs(H_eff(:,13,1)'*v_2) + ...
square_abs(H_eff(:,13,1)'*v_3) + square_abs(H_eff(:,13,2)'*v_4) + square_abs(H_eff(:,13,2)'*v_5) + ...
square_abs(H_eff(:,13,2)'*v_6) + square_abs(H_eff(:,13,3)'*v_8) + square_abs(H_eff(:,13,3)'*v_9) + P_n_eff(13);

T_14 = square_abs(H_eff(:,14,3)'*v_7) + square_abs(H_eff(:,14,1)'*v_1) + square_abs(H_eff(:,14,1)'*v_2) + ...
square_abs(H_eff(:,14,1)'*v_3) + square_abs(H_eff(:,14,2)'*v_4) + square_abs(H_eff(:,14,2)'*v_5) + ...
square_abs(H_eff(:,14,2)'*v_6) + square_abs(H_eff(:,14,3)'*v_8) + square_abs(H_eff(:,14,3)'*v_9) + P_n_eff(14);

T_15 = square_abs(H_eff(:,15,3)'*v_8) + square_abs(H_eff(:,15,1)'*v_1) + square_abs(H_eff(:,15,1)'*v_2) + ...
square_abs(H_eff(:,15,1)'*v_3) + square_abs(H_eff(:,15,2)'*v_4) + square_abs(H_eff(:,15,2)'*v_5) + ...
square_abs(H_eff(:,15,2)'*v_6) + square_abs(H_eff(:,15,3)'*v_7) + square_abs(H_eff(:,15,3)'*v_9) + P_n_eff(15);

T_16 = square_abs(H_eff(:,16,3)'*v_8) + square_abs(H_eff(:,16,1)'*v_1) + square_abs(H_eff(:,16,1)'*v_2) + ...
square_abs(H_eff(:,16,1)'*v_3) + square_abs(H_eff(:,16,2)'*v_4) + square_abs(H_eff(:,16,2)'*v_5) + ...
square_abs(H_eff(:,16,2)'*v_6) + square_abs(H_eff(:,16,3)'*v_7) + square_abs(H_eff(:,16,3)'*v_9) + P_n_eff(16);

T_17 = square_abs(H_eff(:,17,3)'*v_9) + square_abs(H_eff(:,17,1)'*v_1) + square_abs(H_eff(:,17,1)'*v_2) + ...
square_abs(H_eff(:,17,1)'*v_3) + square_abs(H_eff(:,17,2)'*v_4) + square_abs(H_eff(:,17,2)'*v_5) + ...
square_abs(H_eff(:,17,2)'*v_6) + square_abs(H_eff(:,17,3)'*v_7) + square_abs(H_eff(:,17,3)'*v_8) + P_n_eff(17);

T_18 = square_abs(H_eff(:,18,3)'*v_9) + square_abs(H_eff(:,18,1)'*v_1) + square_abs(H_eff(:,18,1)'*v_2) + ...
square_abs(H_eff(:,18,1)'*v_3) + square_abs(H_eff(:,18,2)'*v_4) + square_abs(H_eff(:,18,2)'*v_5) + ...
square_abs(H_eff(:,18,2)'*v_6) + square_abs(H_eff(:,18,3)'*v_7) + square_abs(H_eff(:,18,3)'*v_8) + P_n_eff(18);

%MSE
%E_k Calculation
E_1 = abs(g_k(1))^2*T_1 - 2*real(g_k(1)*H_eff(:,1,1)'*v_1) + 1;
E_2 = abs(g_k(2))^2*T_2 - 2*real(g_k(2)*H_eff(:,2,1)'*v_1) + 1;
E_3 = abs(g_k(3))^2*T_3 - 2*real(g_k(3)*H_eff(:,3,1)'*v_2) + 1;
E_4 = abs(g_k(4))^2*T_4 - 2*real(g_k(4)*H_eff(:,4,1)'*v_2) + 1;
E_5 = abs(g_k(5))^2*T_5 - 2*real(g_k(5)*H_eff(:,5,1)'*v_3) + 1;
E_6 = abs(g_k(6))^2*T_6 - 2*real(g_k(6)*H_eff(:,6,1)'*v_3) + 1;
E_7 = abs(g_k(7))^2*T_7 - 2*real(g_k(7)*H_eff(:,7,2)'*v_4) + 1;
E_8 = abs(g_k(8))^2*T_8 - 2*real(g_k(8)*H_eff(:,8,2)'*v_4) + 1;
E_9 = abs(g_k(9))^2*T_9 - 2*real(g_k(9)*H_eff(:,9,2)'*v_5) + 1;
E_10 = abs(g_k(10))^2*T_10 - 2*real(g_k(10)*H_eff(:,10,2)'*v_5) + 1;
E_11 = abs(g_k(11))^2*T_11 - 2*real(g_k(11)*H_eff(:,11,2)'*v_6) + 1;
E_12 = abs(g_k(12))^2*T_12 - 2*real(g_k(12)*H_eff(:,12,2)'*v_6) + 1;
E_13 = abs(g_k(13))^2*T_13 - 2*real(g_k(13)*H_eff(:,13,3)'*v_7) + 1;
E_14 = abs(g_k(14))^2*T_14 - 2*real(g_k(14)*H_eff(:,14,3)'*v_7) + 1;
E_15 = abs(g_k(15))^2*T_15 - 2*real(g_k(15)*H_eff(:,15,3)'*v_8) + 1;
E_16 = abs(g_k(16))^2*T_16 - 2*real(g_k(16)*H_eff(:,16,3)'*v_8) + 1;
E_17 = abs(g_k(17))^2*T_17 - 2*real(g_k(17)*H_eff(:,17,3)'*v_9) + 1;
E_18 = abs(g_k(18))^2*T_18 - 2*real(g_k(18)*H_eff(:,18,3)'*v_9) + 1;

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

%OBP Covariance Matrix
C1 = R(:,:,1)*R(:,:,1)';
C2 = R(:,:,2)*R(:,:,2)';
C3 = R(:,:,3)*R(:,:,3)';

%Objective Function
object_func = r_g;

%Optimisation
maximize(object_func)

%Constraints    
constraints(1) = r_1 - r_g;
constraints(2) = r_2 - r_g;
constraints(3) = r_3 - r_g;
constraints(4) = r_4 - r_g;
constraints(5) = r_5 - r_g;
constraints(6) = r_6 - r_g;
constraints(7) = r_7 - r_g;
constraints(8) = r_8 - r_g;
constraints(9) = r_9 - r_g;
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

constraints(28) = PAC - (v_1'*G_eff(:,:,1,1)'*D1*G_eff(:,:,1,1)*v_1 + v_2'*G_eff(:,:,1,1)'*D1*G_eff(:,:,1,1)*v_2 + v_3'*G_eff(:,:,1,1)'*D1*G_eff(:,:,1,1)*v_3 + ... 
v_4'*G_eff(:,:,1,2)'*D1*G_eff(:,:,1,2)*v_4 + v_5'*G_eff(:,:,1,2)'*D1*G_eff(:,:,1,2)*v_5 + v_6'*G_eff(:,:,1,2)'*D1*G_eff(:,:,1,2)*v_6 + ...
v_7'*G_eff(:,:,1,3)'*D1*G_eff(:,:,1,3)*v_7 + v_8'*G_eff(:,:,1,3)'*D1*G_eff(:,:,1,3)*v_8 + v_9'*G_eff(:,:,1,3)'*D1*G_eff(:,:,1,3)*v_9 + C1(1,1));

constraints(29) = PAC - (v_1'*G_eff(:,:,1,1)'*D2*G_eff(:,:,1,1)*v_1 + v_2'*G_eff(:,:,1,1)'*D2*G_eff(:,:,1,1)*v_2 + v_3'*G_eff(:,:,1,1)'*D2*G_eff(:,:,1,1)*v_3 + ... 
v_4'*G_eff(:,:,1,2)'*D2*G_eff(:,:,1,2)*v_4 + v_5'*G_eff(:,:,1,2)'*D2*G_eff(:,:,1,2)*v_5 + v_6'*G_eff(:,:,1,2)'*D2*G_eff(:,:,1,2)*v_6 + ...
v_7'*G_eff(:,:,1,3)'*D2*G_eff(:,:,1,3)*v_7 + v_8'*G_eff(:,:,1,3)'*D2*G_eff(:,:,1,3)*v_8 + v_9'*G_eff(:,:,1,3)'*D2*G_eff(:,:,1,3)*v_9 + C1(2,2));

constraints(30) = PAC - (v_1'*G_eff(:,:,1,1)'*D3*G_eff(:,:,1,1)*v_1 + v_2'*G_eff(:,:,1,1)'*D3*G_eff(:,:,1,1)*v_2 + v_3'*G_eff(:,:,1,1)'*D3*G_eff(:,:,1,1)*v_3 + ... 
v_4'*G_eff(:,:,1,2)'*D3*G_eff(:,:,1,2)*v_4 + v_5'*G_eff(:,:,1,2)'*D3*G_eff(:,:,1,2)*v_5 + v_6'*G_eff(:,:,1,2)'*D3*G_eff(:,:,1,2)*v_6 + ...
v_7'*G_eff(:,:,1,3)'*D3*G_eff(:,:,1,3)*v_7 + v_8'*G_eff(:,:,1,3)'*D3*G_eff(:,:,1,3)*v_8 + v_9'*G_eff(:,:,1,3)'*D3*G_eff(:,:,1,3)*v_9 + C1(3,3));

constraints(31) = PAC - (v_1'*G_eff(:,:,2,1)'*D1*G_eff(:,:,2,1)*v_1 + v_2'*G_eff(:,:,2,1)'*D1*G_eff(:,:,2,1)*v_2 + v_3'*G_eff(:,:,2,1)'*D1*G_eff(:,:,2,1)*v_3 + ... 
v_4'*G_eff(:,:,2,2)'*D1*G_eff(:,:,2,2)*v_4 + v_5'*G_eff(:,:,2,2)'*D1*G_eff(:,:,2,2)*v_5 + v_6'*G_eff(:,:,2,2)'*D1*G_eff(:,:,2,2)*v_6 + ...
v_7'*G_eff(:,:,2,3)'*D1*G_eff(:,:,2,3)*v_7 + v_8'*G_eff(:,:,2,3)'*D1*G_eff(:,:,2,3)*v_8 + v_9'*G_eff(:,:,2,3)'*D1*G_eff(:,:,2,3)*v_9 + C2(1,1));

constraints(32) = PAC - (v_1'*G_eff(:,:,2,1)'*D2*G_eff(:,:,2,1)*v_1 + v_2'*G_eff(:,:,2,1)'*D2*G_eff(:,:,2,1)*v_2 + v_3'*G_eff(:,:,2,1)'*D2*G_eff(:,:,2,1)*v_3 + ... 
v_4'*G_eff(:,:,2,2)'*D2*G_eff(:,:,2,2)*v_4 + v_5'*G_eff(:,:,2,2)'*D2*G_eff(:,:,2,2)*v_5 + v_6'*G_eff(:,:,2,2)'*D2*G_eff(:,:,2,2)*v_6 + ...
v_7'*G_eff(:,:,2,3)'*D2*G_eff(:,:,2,3)*v_7 + v_8'*G_eff(:,:,2,3)'*D2*G_eff(:,:,2,3)*v_8 + v_9'*G_eff(:,:,2,3)'*D2*G_eff(:,:,2,3)*v_9 + C2(2,2));

constraints(33) = PAC - (v_1'*G_eff(:,:,2,1)'*D3*G_eff(:,:,2,1)*v_1 + v_2'*G_eff(:,:,2,1)'*D3*G_eff(:,:,2,1)*v_2 + v_3'*G_eff(:,:,2,1)'*D3*G_eff(:,:,2,1)*v_3 + ... 
v_4'*G_eff(:,:,2,2)'*D3*G_eff(:,:,2,2)*v_4 + v_5'*G_eff(:,:,2,2)'*D3*G_eff(:,:,2,2)*v_5 + v_6'*G_eff(:,:,2,2)'*D3*G_eff(:,:,2,2)*v_6 + ...
v_7'*G_eff(:,:,2,3)'*D3*G_eff(:,:,2,3)*v_7 + v_8'*G_eff(:,:,2,3)'*D3*G_eff(:,:,2,3)*v_8 + v_9'*G_eff(:,:,2,3)'*D3*G_eff(:,:,2,3)*v_9 + C2(3,3));

constraints(34) = PAC - (v_1'*G_eff(:,:,3,1)'*D1*G_eff(:,:,3,1)*v_1 + v_2'*G_eff(:,:,3,1)'*D1*G_eff(:,:,3,1)*v_2 + v_3'*G_eff(:,:,3,1)'*D1*G_eff(:,:,3,1)*v_3 + ... 
v_4'*G_eff(:,:,3,2)'*D1*G_eff(:,:,3,2)*v_4 + v_5'*G_eff(:,:,3,2)'*D1*G_eff(:,:,3,2)*v_5 + v_6'*G_eff(:,:,3,2)'*D1*G_eff(:,:,3,2)*v_6 + ...
v_7'*G_eff(:,:,3,3)'*D1*G_eff(:,:,3,3)*v_7 + v_8'*G_eff(:,:,3,3)'*D1*G_eff(:,:,3,3)*v_8 + v_9'*G_eff(:,:,3,3)'*D1*G_eff(:,:,3,3)*v_9 + C3(1,1));

constraints(35) = PAC - (v_1'*G_eff(:,:,3,1)'*D2*G_eff(:,:,3,1)*v_1 + v_2'*G_eff(:,:,3,1)'*D2*G_eff(:,:,3,1)*v_2 + v_3'*G_eff(:,:,3,1)'*D2*G_eff(:,:,3,1)*v_3 + ... 
v_4'*G_eff(:,:,3,2)'*D2*G_eff(:,:,3,2)*v_4 + v_5'*G_eff(:,:,3,2)'*D2*G_eff(:,:,3,2)*v_5 + v_6'*G_eff(:,:,3,2)'*D2*G_eff(:,:,3,2)*v_6 + ...
v_7'*G_eff(:,:,3,3)'*D2*G_eff(:,:,3,3)*v_7 + v_8'*G_eff(:,:,3,3)'*D2*G_eff(:,:,3,3)*v_8 + v_9'*G_eff(:,:,3,3)'*D2*G_eff(:,:,3,3)*v_9 + C3(2,2));

constraints(36) = PAC - (v_1'*G_eff(:,:,3,1)'*D3*G_eff(:,:,3,1)*v_1 + v_2'*G_eff(:,:,3,1)'*D3*G_eff(:,:,3,1)*v_2 + v_3'*G_eff(:,:,3,1)'*D3*G_eff(:,:,3,1)*v_3 + ... 
v_4'*G_eff(:,:,3,2)'*D3*G_eff(:,:,3,2)*v_4 + v_5'*G_eff(:,:,3,2)'*D3*G_eff(:,:,3,2)*v_5 + v_6'*G_eff(:,:,3,2)'*D3*G_eff(:,:,3,2)*v_6 + ...
v_7'*G_eff(:,:,3,3)'*D3*G_eff(:,:,3,3)*v_7 + v_8'*G_eff(:,:,3,3)'*D3*G_eff(:,:,3,3)*v_8 + v_9'*G_eff(:,:,3,3)'*D3*G_eff(:,:,3,3)*v_9 + C3(3,3));


subject to
    constraints >= zeros(size(constraints))

cvx_end

%Result
rate = object_func;
