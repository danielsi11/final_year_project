function [rate,p_sc,p_1,p_2,p_3,p_4,p_5,p_6,p_7,p_8,p_9] = RS_perfect_optimisation(H,PAC,g_k,g_sck,W_k,W_sck)
    
[N,K] = size(H); 

D1 = zeros(N,N);
D2 = zeros(N,N);
D3 = zeros(N,N);
D4 = zeros(N,N);
D5 = zeros(N,N);
D6 = zeros(N,N);
D7 = zeros(N,N);
D8 = zeros(N,N);
D9 = zeros(N,N);

D1(1,1) = 1;
D2(2,2) = 1;
D3(3,3) = 1;
D4(4,4) = 1;
D5(5,5) = 1;
D6(6,6) = 1;
D7(7,7) = 1;
D8(8,8) = 1;
D9(9,9) = 1;

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
variable p_sc(N,1) complex
variable SC_1 
variable SC_2 
variable SC_3
variable SC_4 
variable SC_5 
variable SC_6
variable SC_7
variable SC_8
variable SC_9
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
T_1 = square_abs(H(:,1)'*p_1) + square_abs(H(:,1)'*p_2) + square_abs(H(:,1)'*p_3) + ... 
square_abs(H(:,1)'*p_4) + square_abs(H(:,1)'*p_5) + square_abs(H(:,1)'*p_6) + ...
square_abs(H(:,1)'*p_7) + square_abs(H(:,1)'*p_8) + square_abs(H(:,1)'*p_9) + 1;

T_2 = square_abs(H(:,2)'*p_1) + square_abs(H(:,2)'*p_2) + square_abs(H(:,2)'*p_3) + ...
square_abs(H(:,2)'*p_4) + square_abs(H(:,2)'*p_5) + square_abs(H(:,2)'*p_6) + ...
square_abs(H(:,2)'*p_7) + square_abs(H(:,2)'*p_8) + square_abs(H(:,2)'*p_9) + 1;

T_3 = square_abs(H(:,3)'*p_2) + square_abs(H(:,3)'*p_1) + square_abs(H(:,3)'*p_3) + ...
square_abs(H(:,3)'*p_4) + square_abs(H(:,3)'*p_5) + square_abs(H(:,3)'*p_6) + ...
square_abs(H(:,3)'*p_7) + square_abs(H(:,3)'*p_8) + square_abs(H(:,3)'*p_9) + 1;

T_4 = square_abs(H(:,4)'*p_2) + square_abs(H(:,4)'*p_1) + square_abs(H(:,4)'*p_3) + ...
square_abs(H(:,4)'*p_4) + square_abs(H(:,4)'*p_5) + square_abs(H(:,4)'*p_6) + ...
square_abs(H(:,4)'*p_7) + square_abs(H(:,4)'*p_8) + square_abs(H(:,4)'*p_9) + 1;

T_5 = square_abs(H(:,5)'*p_3) + square_abs(H(:,5)'*p_1) + square_abs(H(:,5)'*p_2) + ...
square_abs(H(:,5)'*p_4) + square_abs(H(:,5)'*p_5) + square_abs(H(:,5)'*p_6) + ...
square_abs(H(:,5)'*p_7) + square_abs(H(:,5)'*p_8) + square_abs(H(:,5)'*p_9) + 1;

T_6 = square_abs(H(:,6)'*p_3) + square_abs(H(:,6)'*p_1) + square_abs(H(:,6)'*p_2) + ...
square_abs(H(:,6)'*p_4) + square_abs(H(:,6)'*p_5) + square_abs(H(:,6)'*p_6) + ...
square_abs(H(:,6)'*p_7) + square_abs(H(:,6)'*p_8) + square_abs(H(:,6)'*p_9) + 1;

T_7 = square_abs(H(:,7)'*p_4) + square_abs(H(:,7)'*p_1) + square_abs(H(:,7)'*p_2) + ...
square_abs(H(:,7)'*p_3) + square_abs(H(:,7)'*p_5) + square_abs(H(:,7)'*p_6) + ...
square_abs(H(:,7)'*p_7) + square_abs(H(:,7)'*p_8) + square_abs(H(:,7)'*p_9) + 1;

T_8 = square_abs(H(:,8)'*p_4) + square_abs(H(:,8)'*p_1) + square_abs(H(:,8)'*p_2) + ...
square_abs(H(:,8)'*p_3) + square_abs(H(:,8)'*p_5) + square_abs(H(:,8)'*p_6) + ...
square_abs(H(:,8)'*p_7) + square_abs(H(:,8)'*p_8) + square_abs(H(:,8)'*p_9) + 1;

T_9 = square_abs(H(:,9)'*p_5) + square_abs(H(:,9)'*p_1) + square_abs(H(:,9)'*p_2) + ...
square_abs(H(:,9)'*p_3) + square_abs(H(:,9)'*p_4) + square_abs(H(:,9)'*p_6) + ...
square_abs(H(:,9)'*p_7) + square_abs(H(:,9)'*p_8) + square_abs(H(:,9)'*p_9) + 1;

T_10 = square_abs(H(:,10)'*p_5) + square_abs(H(:,10)'*p_1) + square_abs(H(:,10)'*p_2) + ...
square_abs(H(:,10)'*p_3) + square_abs(H(:,10)'*p_4) + square_abs(H(:,10)'*p_6) + ...
square_abs(H(:,10)'*p_7) + square_abs(H(:,10)'*p_8) + square_abs(H(:,10)'*p_9) + 1;

T_11 = square_abs(H(:,11)'*p_6) + square_abs(H(:,11)'*p_1) + square_abs(H(:,11)'*p_2) + ...
square_abs(H(:,11)'*p_3) + square_abs(H(:,11)'*p_4) + square_abs(H(:,11)'*p_5) + ...
square_abs(H(:,11)'*p_7) + square_abs(H(:,11)'*p_8) + square_abs(H(:,11)'*p_9) + 1;

T_12 = square_abs(H(:,12)'*p_6) + square_abs(H(:,12)'*p_1) + square_abs(H(:,12)'*p_2) + ...
square_abs(H(:,12)'*p_3) + square_abs(H(:,12)'*p_4) + square_abs(H(:,12)'*p_5) + ...
square_abs(H(:,12)'*p_7) + square_abs(H(:,12)'*p_8) + square_abs(H(:,12)'*p_9) + 1;

T_13 = square_abs(H(:,13)'*p_7) + square_abs(H(:,13)'*p_1) + square_abs(H(:,13)'*p_2) + ...
square_abs(H(:,13)'*p_3) + square_abs(H(:,13)'*p_4) + square_abs(H(:,13)'*p_5) + ...
square_abs(H(:,13)'*p_6) + square_abs(H(:,13)'*p_8) + square_abs(H(:,13)'*p_9) + 1;

T_14 = square_abs(H(:,14)'*p_7) + square_abs(H(:,14)'*p_1) + square_abs(H(:,14)'*p_2) + ...
square_abs(H(:,14)'*p_3) + square_abs(H(:,14)'*p_4) + square_abs(H(:,14)'*p_5) + ...
square_abs(H(:,14)'*p_6) + square_abs(H(:,14)'*p_8) + square_abs(H(:,14)'*p_9) + 1;

T_15 = square_abs(H(:,15)'*p_8) + square_abs(H(:,15)'*p_1) + square_abs(H(:,15)'*p_2) + ...
square_abs(H(:,15)'*p_3) + square_abs(H(:,15)'*p_4) + square_abs(H(:,15)'*p_5) + ...
square_abs(H(:,15)'*p_6) + square_abs(H(:,15)'*p_7) + square_abs(H(:,15)'*p_9) + 1;

T_16 = square_abs(H(:,16)'*p_8) + square_abs(H(:,16)'*p_1) + square_abs(H(:,16)'*p_2) + ...
square_abs(H(:,16)'*p_3) + square_abs(H(:,16)'*p_4) + square_abs(H(:,16)'*p_5) + ...
square_abs(H(:,16)'*p_6) + square_abs(H(:,16)'*p_7) + square_abs(H(:,16)'*p_9) + 1;

T_17 = square_abs(H(:,17)'*p_9) + square_abs(H(:,17)'*p_1) + square_abs(H(:,17)'*p_2) + ...
square_abs(H(:,17)'*p_3) + square_abs(H(:,17)'*p_4) + square_abs(H(:,17)'*p_5) + ...
square_abs(H(:,17)'*p_6) + square_abs(H(:,17)'*p_7) + square_abs(H(:,17)'*p_8) + 1;

T_18 = square_abs(H(:,18)'*p_9) + square_abs(H(:,18)'*p_1) + square_abs(H(:,18)'*p_2) + ...
square_abs(H(:,18)'*p_3) + square_abs(H(:,18)'*p_4) + square_abs(H(:,18)'*p_5) + ...
square_abs(H(:,18)'*p_6) + square_abs(H(:,18)'*p_7) + square_abs(H(:,18)'*p_8) + 1;

%T_sck Calculation
T_sc1 = square_abs(H(:,1)'*p_sc) + T_1;
T_sc2 = square_abs(H(:,2)'*p_sc) + T_2;
T_sc3 = square_abs(H(:,3)'*p_sc) + T_3;
T_sc4 = square_abs(H(:,4)'*p_sc) + T_4;
T_sc5 = square_abs(H(:,5)'*p_sc) + T_5;
T_sc6 = square_abs(H(:,6)'*p_sc) + T_6;
T_sc7 = square_abs(H(:,7)'*p_sc) + T_7;
T_sc8 = square_abs(H(:,8)'*p_sc) + T_8;
T_sc9 = square_abs(H(:,9)'*p_sc) + T_9;
T_sc10 = square_abs(H(:,10)'*p_sc) + T_10;
T_sc11 = square_abs(H(:,11)'*p_sc) + T_11;
T_sc12 = square_abs(H(:,12)'*p_sc) + T_12;
T_sc13 = square_abs(H(:,13)'*p_sc) + T_13;
T_sc14 = square_abs(H(:,14)'*p_sc) + T_14;
T_sc15 = square_abs(H(:,15)'*p_sc) + T_15;
T_sc16 = square_abs(H(:,16)'*p_sc) + T_16;
T_sc17 = square_abs(H(:,17)'*p_sc) + T_17;
T_sc18 = square_abs(H(:,18)'*p_sc) + T_18;

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
E_15 = abs(g_k(15))^2*T_15 - 2*real(g_k(15)*H(:,15)'*p_8) + 1;
E_16 = abs(g_k(16))^2*T_16 - 2*real(g_k(16)*H(:,16)'*p_8) + 1;
E_17 = abs(g_k(17))^2*T_17 - 2*real(g_k(17)*H(:,17)'*p_9) + 1;
E_18 = abs(g_k(18))^2*T_18 - 2*real(g_k(18)*H(:,18)'*p_9) + 1;

%E_sck Calculation
E_sc1 = abs(g_sck(1))^2*T_sc1 - 2*real(g_sck(1)*H(:,1)'*p_sc) + 1;
E_sc2 = abs(g_sck(2))^2*T_sc2 - 2*real(g_sck(2)*H(:,2)'*p_sc) + 1;
E_sc3 = abs(g_sck(3))^2*T_sc3 - 2*real(g_sck(3)*H(:,3)'*p_sc) + 1;
E_sc4 = abs(g_sck(4))^2*T_sc4 - 2*real(g_sck(4)*H(:,4)'*p_sc) + 1;
E_sc5 = abs(g_sck(5))^2*T_sc5 - 2*real(g_sck(5)*H(:,5)'*p_sc) + 1;
E_sc6 = abs(g_sck(6))^2*T_sc6 - 2*real(g_sck(6)*H(:,6)'*p_sc) + 1;
E_sc7 = abs(g_sck(7))^2*T_sc7 - 2*real(g_sck(7)*H(:,7)'*p_sc) + 1;
E_sc8 = abs(g_sck(8))^2*T_sc8 - 2*real(g_sck(8)*H(:,8)'*p_sc) + 1;
E_sc9 = abs(g_sck(9))^2*T_sc9 - 2*real(g_sck(9)*H(:,9)'*p_sc) + 1;
E_sc10 = abs(g_sck(10))^2*T_sc10 - 2*real(g_sck(10)*H(:,10)'*p_sc) + 1;
E_sc11 = abs(g_sck(11))^2*T_sc11 - 2*real(g_sck(11)*H(:,11)'*p_sc) + 1;
E_sc12 = abs(g_sck(12))^2*T_sc12 - 2*real(g_sck(12)*H(:,12)'*p_sc) + 1;
E_sc13 = abs(g_sck(13))^2*T_sc13 - 2*real(g_sck(13)*H(:,13)'*p_sc) + 1;
E_sc14 = abs(g_sck(14))^2*T_sc14 - 2*real(g_sck(14)*H(:,14)'*p_sc) + 1;
E_sc15 = abs(g_sck(15))^2*T_sc15 - 2*real(g_sck(15)*H(:,15)'*p_sc) + 1;
E_sc16 = abs(g_sck(16))^2*T_sc16 - 2*real(g_sck(16)*H(:,16)'*p_sc) + 1;
E_sc17 = abs(g_sck(17))^2*T_sc17 - 2*real(g_sck(17)*H(:,17)'*p_sc) + 1;
E_sc18 = abs(g_sck(18))^2*T_sc18 - 2*real(g_sck(18)*H(:,18)'*p_sc) + 1;

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

%Global Common Stream
X_sc1 = W_sck(1)*E_sc1 - log2(W_sck(1));
X_sc2 = W_sck(2)*E_sc2 - log2(W_sck(2));
X_sc3 = W_sck(3)*E_sc3 - log2(W_sck(3));
X_sc4 = W_sck(4)*E_sc4 - log2(W_sck(4));
X_sc5 = W_sck(5)*E_sc5 - log2(W_sck(5));
X_sc6 = W_sck(6)*E_sc6 - log2(W_sck(6));
X_sc7 = W_sck(7)*E_sc7 - log2(W_sck(7));
X_sc8 = W_sck(8)*E_sc8 - log2(W_sck(8));
X_sc9 = W_sck(9)*E_sc9 - log2(W_sck(9));
X_sc10 = W_sck(10)*E_sc10 - log2(W_sck(10));
X_sc11 = W_sck(11)*E_sc11 - log2(W_sck(11));
X_sc12 = W_sck(12)*E_sc12 - log2(W_sck(12));
X_sc13 = W_sck(13)*E_sc13 - log2(W_sck(13));
X_sc14 = W_sck(14)*E_sc14 - log2(W_sck(14));
X_sc15 = W_sck(15)*E_sc15 - log2(W_sck(15));
X_sc16 = W_sck(16)*E_sc16 - log2(W_sck(16));
X_sc17 = W_sck(17)*E_sc17 - log2(W_sck(17));
X_sc18 = W_sck(18)*E_sc18 - log2(W_sck(18));

%Objective Function
object_func = r_g;

%Optimisation
maximize(object_func)

%Constraints    
constraints(1) = SC_1 + r_1 - r_g;
constraints(2) = SC_2 + r_2 - r_g;
constraints(3) = SC_3 + r_3 - r_g;
constraints(4) = SC_4 + r_4 - r_g;
constraints(5) = SC_5 + r_5 - r_g;
constraints(6) = SC_6 + r_6 - r_g;
constraints(7) = SC_7 + r_7 - r_g;
constraints(8) = SC_8 + r_8 - r_g;
constraints(9) = SC_9 + r_9 - r_g;
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
constraints(28) = 1 - X_sc1 - SC_1 - SC_2 - SC_3 - SC_4 - SC_5 - SC_6 - SC_7 - SC_8 - SC_9;
constraints(29) = 1 - X_sc2 - SC_1 - SC_2 - SC_3 - SC_4 - SC_5 - SC_6 - SC_7 - SC_8 - SC_9;
constraints(30) = 1 - X_sc3 - SC_1 - SC_2 - SC_3 - SC_4 - SC_5 - SC_6 - SC_7 - SC_8 - SC_9;
constraints(31) = 1 - X_sc4 - SC_1 - SC_2 - SC_3 - SC_4 - SC_5 - SC_6 - SC_7 - SC_8 - SC_9;
constraints(32) = 1 - X_sc5 - SC_1 - SC_2 - SC_3 - SC_4 - SC_5 - SC_6 - SC_7 - SC_8 - SC_9;
constraints(33) = 1 - X_sc6 - SC_1 - SC_2 - SC_3 - SC_4 - SC_5 - SC_6 - SC_7 - SC_8 - SC_9;
constraints(34) = 1 - X_sc7 - SC_1 - SC_2 - SC_3 - SC_4 - SC_5 - SC_6 - SC_7 - SC_8 - SC_9;
constraints(35) = 1 - X_sc8 - SC_1 - SC_2 - SC_3 - SC_4 - SC_5 - SC_6 - SC_7 - SC_8 - SC_9;
constraints(36) = 1 - X_sc9 - SC_1 - SC_2 - SC_3 - SC_4 - SC_5 - SC_6 - SC_7 - SC_8 - SC_9;
constraints(37) = 1 - X_sc10 - SC_1 - SC_2 - SC_3 - SC_4 - SC_5 - SC_6 - SC_7 - SC_8 - SC_9;
constraints(38) = 1 - X_sc11 - SC_1 - SC_2 - SC_3 - SC_4 - SC_5 - SC_6 - SC_7 - SC_8 - SC_9;
constraints(39) = 1 - X_sc12 - SC_1 - SC_2 - SC_3 - SC_4 - SC_5 - SC_6 - SC_7 - SC_8 - SC_9;
constraints(40) = 1 - X_sc13 - SC_1 - SC_2 - SC_3 - SC_4 - SC_5 - SC_6 - SC_7 - SC_8 - SC_9;
constraints(41) = 1 - X_sc14 - SC_1 - SC_2 - SC_3 - SC_4 - SC_5 - SC_6 - SC_7 - SC_8 - SC_9;
constraints(42) = 1 - X_sc15 - SC_1 - SC_2 - SC_3 - SC_4 - SC_5 - SC_6 - SC_7 - SC_8 - SC_9;
constraints(43) = 1 - X_sc16 - SC_1 - SC_2 - SC_3 - SC_4 - SC_5 - SC_6 - SC_7 - SC_8 - SC_9;
constraints(44) = 1 - X_sc17 - SC_1 - SC_2 - SC_3 - SC_4 - SC_5 - SC_6 - SC_7 - SC_8 - SC_9;
constraints(45) = 1 - X_sc18 - SC_1 - SC_2 - SC_3 - SC_4 - SC_5 - SC_6 - SC_7 - SC_8 - SC_9;
constraints(46) = SC_1;
constraints(47) = SC_2;
constraints(48) = SC_3;
constraints(49) = SC_4;
constraints(50) = SC_5;
constraints(51) = SC_6;
constraints(52) = SC_7;
constraints(53) = SC_8;
constraints(54) = SC_9;
constraints(55) = PAC - (p_sc'*D1*p_sc + p_1'*D1*p_1 + p_2'*D1*p_2 + p_3'*D1*p_3 +...
p_4'*D1*p_4 + p_5'*D1*p_5 + p_6'*D1*p_6 + p_7'*D1*p_7 + p_8'*D1*p_8 + p_9'*D1*p_9);
constraints(56) = PAC - (p_sc'*D2*p_sc + p_1'*D2*p_1 + p_2'*D2*p_2 + p_3'*D2*p_3 +...
p_4'*D2*p_4 + p_5'*D2*p_5 + p_6'*D2*p_6 + p_7'*D2*p_7 + p_8'*D2*p_8 + p_9'*D2*p_9);
constraints(57) = PAC - (p_sc'*D3*p_sc + p_1'*D3*p_1 + p_2'*D3*p_2 + p_3'*D3*p_3 +...
p_4'*D3*p_4 + p_5'*D3*p_5 + p_6'*D3*p_6 + p_7'*D3*p_7 + p_8'*D3*p_8 + p_9'*D3*p_9);
constraints(58) = PAC - (p_sc'*D4*p_sc + p_1'*D4*p_1 + p_2'*D4*p_2 + p_3'*D4*p_3 +...
p_4'*D4*p_4 + p_5'*D4*p_5 + p_6'*D4*p_6 + p_7'*D4*p_7 + p_8'*D4*p_8 + p_9'*D4*p_9);
constraints(59) = PAC - (p_sc'*D5*p_sc + p_1'*D5*p_1 + p_2'*D5*p_2 + p_3'*D5*p_3 +...
p_4'*D5*p_4 + p_5'*D5*p_5 + p_6'*D5*p_6 + p_7'*D5*p_7 + p_8'*D5*p_8 + p_9'*D5*p_9);
constraints(60) = PAC - (p_sc'*D6*p_sc + p_1'*D6*p_1 + p_2'*D6*p_2 + p_3'*D6*p_3 +...
p_4'*D6*p_4 + p_5'*D6*p_5 + p_6'*D6*p_6 + p_7'*D6*p_7 + p_8'*D6*p_8 + p_9'*D6*p_9);
constraints(61) = PAC - (p_sc'*D7*p_sc + p_1'*D7*p_1 + p_2'*D7*p_2 + p_3'*D7*p_3 +...
p_4'*D7*p_4 + p_5'*D7*p_5 + p_6'*D7*p_6 + p_7'*D7*p_7 + p_8'*D7*p_8 + p_9'*D7*p_9);
constraints(62) = PAC - (p_sc'*D8*p_sc + p_1'*D8*p_1 + p_2'*D8*p_2 + p_3'*D8*p_3 +...
p_4'*D8*p_4 + p_5'*D8*p_5 + p_6'*D8*p_6 + p_7'*D8*p_7 + p_8'*D8*p_8 + p_9'*D8*p_9);
constraints(63) = PAC - (p_sc'*D9*p_sc + p_1'*D9*p_1 + p_2'*D9*p_2 + p_3'*D9*p_3 +...
p_4'*D9*p_4 + p_5'*D9*p_5 + p_6'*D9*p_6 + p_7'*D9*p_7 + p_8'*D9*p_8 + p_9'*D9*p_9);

subject to
    constraints >= zeros(size(constraints))

cvx_end

%Result
rate = object_func;
