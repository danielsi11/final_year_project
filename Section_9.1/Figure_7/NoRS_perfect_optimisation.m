function [rate,p_1,p_2,p_3,p_4,p_5,p_6,p_7] = NoRS_perfect_optimisation(H,PAC,g_k,W_k)
    
[N,K] = size(H); 

D1 = zeros(N,N);
D2 = zeros(N,N);
D3 = zeros(N,N);
D4 = zeros(N,N);
D5 = zeros(N,N);
D6 = zeros(N,N);
D7 = zeros(N,N);

D1(1,1) = 1;
D2(2,2) = 1;
D3(3,3) = 1;
D4(4,4) = 1;
D5(5,5) = 1;
D6(6,6) = 1;
D7(7,7) = 1;

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
variable r_1 
variable r_2 
variable r_3
variable r_4 
variable r_5 
variable r_6
variable r_7
variable r_g

expression constraints(1,28);

%Received Power
%T_k Calculation
T_1 = square_abs(H(:,1)'*p_1) + square_abs(H(:,1)'*p_2) + square_abs(H(:,1)'*p_3) + ...
square_abs(H(:,1)'*p_4) + square_abs(H(:,1)'*p_5) + square_abs(H(:,1)'*p_6) + square_abs(H(:,1)'*p_7) + 1;

T_2 = square_abs(H(:,2)'*p_1) + square_abs(H(:,2)'*p_2) + square_abs(H(:,2)'*p_3) + ...
square_abs(H(:,2)'*p_4) + square_abs(H(:,2)'*p_5) + square_abs(H(:,2)'*p_6) + square_abs(H(:,2)'*p_7) + 1;

T_3 = square_abs(H(:,3)'*p_2) + square_abs(H(:,3)'*p_1) + square_abs(H(:,3)'*p_3) + ...
square_abs(H(:,3)'*p_4) + square_abs(H(:,3)'*p_5) + square_abs(H(:,3)'*p_6) + square_abs(H(:,3)'*p_7) + 1;

T_4 = square_abs(H(:,4)'*p_2) + square_abs(H(:,4)'*p_1) + square_abs(H(:,4)'*p_3) + ...
square_abs(H(:,4)'*p_4) + square_abs(H(:,4)'*p_5) + square_abs(H(:,4)'*p_6) + square_abs(H(:,4)'*p_7) + 1;

T_5 = square_abs(H(:,5)'*p_3) + square_abs(H(:,5)'*p_1) + square_abs(H(:,5)'*p_2) + ...
square_abs(H(:,5)'*p_4) + square_abs(H(:,5)'*p_5) + square_abs(H(:,5)'*p_6) + square_abs(H(:,5)'*p_7) + 1;

T_6 = square_abs(H(:,6)'*p_3) + square_abs(H(:,6)'*p_1) + square_abs(H(:,6)'*p_2) + ...
square_abs(H(:,6)'*p_4) + square_abs(H(:,6)'*p_5) + square_abs(H(:,6)'*p_6) + square_abs(H(:,6)'*p_7) + 1;

T_7 = square_abs(H(:,7)'*p_4) + square_abs(H(:,7)'*p_1) + square_abs(H(:,7)'*p_2) + ...
square_abs(H(:,7)'*p_3) + square_abs(H(:,7)'*p_5) + square_abs(H(:,7)'*p_6) + square_abs(H(:,7)'*p_7) + 1;

T_8 = square_abs(H(:,8)'*p_4) + square_abs(H(:,8)'*p_1) + square_abs(H(:,8)'*p_2) + ...
square_abs(H(:,8)'*p_3) + square_abs(H(:,8)'*p_5) + square_abs(H(:,8)'*p_6) + square_abs(H(:,8)'*p_7) + 1;

T_9 = square_abs(H(:,9)'*p_5) + square_abs(H(:,9)'*p_1) + square_abs(H(:,9)'*p_2) + ...
square_abs(H(:,9)'*p_3) + square_abs(H(:,9)'*p_4) + square_abs(H(:,9)'*p_6) + square_abs(H(:,9)'*p_7) + 1;

T_10 = square_abs(H(:,10)'*p_5) + square_abs(H(:,10)'*p_1) + square_abs(H(:,10)'*p_2) + ...
square_abs(H(:,10)'*p_3) + square_abs(H(:,10)'*p_4) + square_abs(H(:,10)'*p_6) + square_abs(H(:,10)'*p_7) + 1;

T_11 = square_abs(H(:,11)'*p_6) + square_abs(H(:,11)'*p_1) + square_abs(H(:,11)'*p_2) + ...
square_abs(H(:,11)'*p_3) + square_abs(H(:,11)'*p_4) + square_abs(H(:,11)'*p_5) + square_abs(H(:,11)'*p_7) + 1;

T_12 = square_abs(H(:,12)'*p_6) + square_abs(H(:,12)'*p_1) + square_abs(H(:,12)'*p_2) + ...
square_abs(H(:,12)'*p_3) + square_abs(H(:,12)'*p_4) + square_abs(H(:,12)'*p_5) + square_abs(H(:,12)'*p_7) + 1;

T_13 = square_abs(H(:,13)'*p_7) + square_abs(H(:,13)'*p_1) + square_abs(H(:,13)'*p_2) + ...
square_abs(H(:,13)'*p_3) + square_abs(H(:,13)'*p_4) + square_abs(H(:,13)'*p_5) + square_abs(H(:,13)'*p_6) + 1;

T_14 = square_abs(H(:,14)'*p_7) + square_abs(H(:,14)'*p_1) + square_abs(H(:,14)'*p_2) + ...
square_abs(H(:,14)'*p_3) + square_abs(H(:,14)'*p_4) + square_abs(H(:,14)'*p_5) + square_abs(H(:,14)'*p_6) + 1;

%MSE
%E_k Calculation
E_1 = abs(g_k(1))^2*T_1 - 2*real(g_k(1)*H(:,1)'*p_1) + 1;
E_2 = abs(g_k(2))^2*T_2 - 2*real(g_k(2)*H(:,2)'*p_1) + 1;
E_3 = abs(g_k(3))^2*T_3 - 2*real(g_k(3)*H(:,3)'*p_2) + 1;
E_4 = abs(g_k(4))^2*T_4 - 2*real(g_k(4)*H(:,4)'*p_2) + 1;
E_5 = abs(g_k(5))^2*T_5 - 2*real(g_k(5)*H(:,5)'*p_3) + 1;
E_6 = abs(g_k(6))^2*T_6 - 2*real(g_k(6)*H(:,6)'*p_3) + 1;
E_7 = abs(g_k(7))^2*T_7 - 2*real(g_k(7)*H(:,7)'*p_4) + 1;
E_8 = abs(g_k(8))^2*T_8 - 2*real(g_k(8)*H(:,8)'*p_4) + 1;
E_9 = abs(g_k(9))^2*T_9 - 2*real(g_k(9)*H(:,9)'*p_5) + 1;
E_10 = abs(g_k(10))^2*T_10 - 2*real(g_k(10)*H(:,10)'*p_5) + 1;
E_11 = abs(g_k(11))^2*T_11 - 2*real(g_k(11)*H(:,11)'*p_6) + 1;
E_12 = abs(g_k(12))^2*T_12 - 2*real(g_k(12)*H(:,12)'*p_6) + 1;
E_13 = abs(g_k(13))^2*T_13 - 2*real(g_k(13)*H(:,13)'*p_7) + 1;
E_14 = abs(g_k(14))^2*T_14 - 2*real(g_k(14)*H(:,14)'*p_7) + 1;

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
constraints(8) = 1 - X_1 - r_1;
constraints(9) = 1 - X_2 - r_1;
constraints(10) = 1 - X_3 - r_2;
constraints(11) = 1 - X_4 - r_2;
constraints(12) = 1 - X_5 - r_3;
constraints(13) = 1 - X_6 - r_3;
constraints(14) = 1 - X_7 - r_4;
constraints(15) = 1 - X_8 - r_4;
constraints(16) = 1 - X_9 - r_5;
constraints(17) = 1 - X_10 - r_5;
constraints(18) = 1 - X_11 - r_6;
constraints(19) = 1 - X_12 - r_6;
constraints(20) = 1 - X_13 - r_7;
constraints(21) = 1 - X_14 - r_7;
constraints(22) = PAC - (p_1'*D1*p_1 + p_2'*D1*p_2 + p_3'*D1*p_3 + ...
p_4'*D1*p_4 + p_5'*D1*p_5 + p_6'*D1*p_6 + p_7'*D1*p_7);
constraints(23) = PAC - (p_1'*D2*p_1 + p_2'*D2*p_2 + p_3'*D2*p_3 + ...
p_4'*D2*p_4 + p_5'*D2*p_5 + p_6'*D2*p_6 + p_7'*D2*p_7);
constraints(24) = PAC - (p_1'*D3*p_1 + p_2'*D3*p_2 + p_3'*D3*p_3 + ...
p_4'*D3*p_4 + p_5'*D3*p_5 + p_6'*D3*p_6 + p_7'*D3*p_7);
constraints(25) = PAC - (p_1'*D4*p_1 + p_2'*D4*p_2 + p_3'*D4*p_3 + ...
p_4'*D4*p_4 + p_5'*D4*p_5 + p_6'*D4*p_6 + p_7'*D4*p_7);
constraints(26) = PAC - (p_1'*D5*p_1 + p_2'*D5*p_2 + p_3'*D5*p_3 + ...
p_4'*D5*p_4 + p_5'*D5*p_5 + p_6'*D5*p_6 + p_7'*D5*p_7);
constraints(27) = PAC - (p_1'*D6*p_1 + p_2'*D6*p_2 + p_3'*D6*p_3 + ...
p_4'*D6*p_4 + p_5'*D6*p_5 + p_6'*D6*p_6 + p_7'*D6*p_7);
constraints(28) = PAC - (p_1'*D7*p_1 + p_2'*D7*p_2 + p_3'*D7*p_3 + ...
p_4'*D7*p_4 + p_5'*D7*p_5 + p_6'*D7*p_6 + p_7'*D7*p_7);

subject to
    constraints >= zeros(size(constraints))

cvx_end

%Result
rate = object_func;
