function [rate,p_1,p_2,p_3,p_4,p_5,p_6,p_7,p_8,p_9] = ...
NoRS_imperfect_optimisation(H_est,PAC,t_k,psi_k,f_k,v_k,u_k,P_n_eff,G)
    
[N,K,L] = size(H_est); 

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

%Rate-WMMSE Relationship
%Private Stream
X_1 = quad_form(p_1,psi_k(:,:,1,1)) + quad_form(p_2,psi_k(:,:,1,1)) + quad_form(p_3,psi_k(:,:,1,1)) + ...
quad_form(p_4,psi_k(:,:,1,2)) + quad_form(p_5,psi_k(:,:,1,2)) + quad_form(p_6,psi_k(:,:,1,2)) + ...
quad_form(p_7,psi_k(:,:,1,3)) + quad_form(p_8,psi_k(:,:,1,3)) + quad_form(p_9,psi_k(:,:,1,3)) + ...
P_n_eff(1)*t_k(1) - 2*real(f_k(:,1)'*p_1) + u_k(1) - v_k(1);

X_2 = quad_form(p_1,psi_k(:,:,2,1)) + quad_form(p_2,psi_k(:,:,2,1)) + quad_form(p_3,psi_k(:,:,2,1)) + ...
quad_form(p_4,psi_k(:,:,2,2)) + quad_form(p_5,psi_k(:,:,2,2)) + quad_form(p_6,psi_k(:,:,2,2)) + ...
quad_form(p_7,psi_k(:,:,2,3)) + quad_form(p_8,psi_k(:,:,2,3)) + quad_form(p_9,psi_k(:,:,2,3)) + ...
P_n_eff(2)*t_k(2) - 2*real(f_k(:,2)'*p_1) + u_k(2) - v_k(2);

X_3 = quad_form(p_1,psi_k(:,:,3,1)) + quad_form(p_2,psi_k(:,:,3,1)) + quad_form(p_3,psi_k(:,:,3,1)) + ...
quad_form(p_4,psi_k(:,:,3,2)) + quad_form(p_5,psi_k(:,:,3,2)) + quad_form(p_6,psi_k(:,:,3,2)) + ...
quad_form(p_7,psi_k(:,:,3,3)) + quad_form(p_8,psi_k(:,:,3,3)) + quad_form(p_9,psi_k(:,:,3,3)) + ...
P_n_eff(3)*t_k(3) - 2*real(f_k(:,3)'*p_2) + u_k(3) - v_k(3);

X_4 = quad_form(p_1,psi_k(:,:,4,1)) + quad_form(p_2,psi_k(:,:,4,1)) + quad_form(p_3,psi_k(:,:,4,1)) + ...
quad_form(p_4,psi_k(:,:,4,2)) + quad_form(p_5,psi_k(:,:,4,2)) + quad_form(p_6,psi_k(:,:,4,2)) + ...
quad_form(p_7,psi_k(:,:,4,3)) + quad_form(p_8,psi_k(:,:,4,3)) + quad_form(p_9,psi_k(:,:,4,3)) + ...
P_n_eff(4)*t_k(4) - 2*real(f_k(:,4)'*p_2) + u_k(4) - v_k(4);

X_5 = quad_form(p_1,psi_k(:,:,5,1)) + quad_form(p_2,psi_k(:,:,5,1)) + quad_form(p_3,psi_k(:,:,5,1)) + ...
quad_form(p_4,psi_k(:,:,5,2)) + quad_form(p_5,psi_k(:,:,5,2)) + quad_form(p_6,psi_k(:,:,5,2)) + ...
quad_form(p_7,psi_k(:,:,5,3)) + quad_form(p_8,psi_k(:,:,5,3)) + quad_form(p_9,psi_k(:,:,5,3)) + ...
P_n_eff(5)*t_k(5) - 2*real(f_k(:,5)'*p_3) + u_k(5) - v_k(5);

X_6 = quad_form(p_1,psi_k(:,:,6,1)) + quad_form(p_2,psi_k(:,:,6,1)) + quad_form(p_3,psi_k(:,:,6,1)) + ...
quad_form(p_4,psi_k(:,:,6,2)) + quad_form(p_5,psi_k(:,:,6,2)) + quad_form(p_6,psi_k(:,:,6,2)) + ...
quad_form(p_7,psi_k(:,:,6,3)) + quad_form(p_8,psi_k(:,:,6,3)) + quad_form(p_9,psi_k(:,:,6,3)) + ...
P_n_eff(6)*t_k(6) - 2*real(f_k(:,6)'*p_3) + u_k(6) - v_k(6);

X_7 = quad_form(p_1,psi_k(:,:,7,1)) + quad_form(p_2,psi_k(:,:,7,1)) + quad_form(p_3,psi_k(:,:,7,1)) + ...
quad_form(p_4,psi_k(:,:,7,2)) + quad_form(p_5,psi_k(:,:,7,2)) + quad_form(p_6,psi_k(:,:,7,2)) + ...
quad_form(p_7,psi_k(:,:,7,3)) + quad_form(p_8,psi_k(:,:,7,3)) + quad_form(p_9,psi_k(:,:,7,3)) + ...
P_n_eff(7)*t_k(7) - 2*real(f_k(:,7)'*p_4) + u_k(7) - v_k(7);

X_8 = quad_form(p_1,psi_k(:,:,8,1)) + quad_form(p_2,psi_k(:,:,8,1)) + quad_form(p_3,psi_k(:,:,8,1)) + ...
quad_form(p_4,psi_k(:,:,8,2)) + quad_form(p_5,psi_k(:,:,8,2)) + quad_form(p_6,psi_k(:,:,8,2)) + ...
quad_form(p_7,psi_k(:,:,8,3)) + quad_form(p_8,psi_k(:,:,8,3)) + quad_form(p_9,psi_k(:,:,8,3)) + ...
P_n_eff(8)*t_k(8) - 2*real(f_k(:,8)'*p_4) + u_k(8) - v_k(8);

X_9 = quad_form(p_1,psi_k(:,:,9,1)) + quad_form(p_2,psi_k(:,:,9,1)) + quad_form(p_3,psi_k(:,:,9,1)) + ...
quad_form(p_4,psi_k(:,:,9,2)) + quad_form(p_5,psi_k(:,:,9,2)) + quad_form(p_6,psi_k(:,:,9,2)) + ...
quad_form(p_7,psi_k(:,:,9,3)) + quad_form(p_8,psi_k(:,:,9,3)) + quad_form(p_9,psi_k(:,:,9,3)) + ...
P_n_eff(9)*t_k(9) - 2*real(f_k(:,9)'*p_5) + u_k(9) - v_k(9);

X_10 = quad_form(p_1,psi_k(:,:,10,1)) + quad_form(p_2,psi_k(:,:,10,1)) + quad_form(p_3,psi_k(:,:,10,1)) + ...
quad_form(p_4,psi_k(:,:,10,2)) + quad_form(p_5,psi_k(:,:,10,2)) + quad_form(p_6,psi_k(:,:,10,2)) + ...
quad_form(p_7,psi_k(:,:,10,3)) + quad_form(p_8,psi_k(:,:,10,3)) + quad_form(p_9,psi_k(:,:,10,3)) + ...
P_n_eff(10)*t_k(10) - 2*real(f_k(:,10)'*p_5) + u_k(10) - v_k(10);

X_11 = quad_form(p_1,psi_k(:,:,11,1)) + quad_form(p_2,psi_k(:,:,11,1)) + quad_form(p_3,psi_k(:,:,11,1)) + ...
quad_form(p_4,psi_k(:,:,11,2)) + quad_form(p_5,psi_k(:,:,11,2)) + quad_form(p_6,psi_k(:,:,11,2)) + ...
quad_form(p_7,psi_k(:,:,11,3)) + quad_form(p_8,psi_k(:,:,11,3)) + quad_form(p_9,psi_k(:,:,11,3)) + ...
P_n_eff(11)*t_k(11) - 2*real(f_k(:,11)'*p_6) + u_k(11) - v_k(11);

X_12 = quad_form(p_1,psi_k(:,:,12,1)) + quad_form(p_2,psi_k(:,:,12,1)) + quad_form(p_3,psi_k(:,:,12,1)) + ...
quad_form(p_4,psi_k(:,:,12,2)) + quad_form(p_5,psi_k(:,:,12,2)) + quad_form(p_6,psi_k(:,:,12,2)) + ...
quad_form(p_7,psi_k(:,:,12,3)) + quad_form(p_8,psi_k(:,:,12,3)) + quad_form(p_9,psi_k(:,:,12,3)) + ...
P_n_eff(12)*t_k(12) - 2*real(f_k(:,12)'*p_6) + u_k(12) - v_k(12);

X_13 = quad_form(p_1,psi_k(:,:,13,1)) + quad_form(p_2,psi_k(:,:,13,1)) + quad_form(p_3,psi_k(:,:,13,1)) + ...
quad_form(p_4,psi_k(:,:,13,2)) + quad_form(p_5,psi_k(:,:,13,2)) + quad_form(p_6,psi_k(:,:,13,2)) + ...
quad_form(p_7,psi_k(:,:,13,3)) + quad_form(p_8,psi_k(:,:,13,3)) + quad_form(p_9,psi_k(:,:,13,3)) + ...
P_n_eff(13)*t_k(13) - 2*real(f_k(:,13)'*p_7) + u_k(13) - v_k(13);

X_14 = quad_form(p_1,psi_k(:,:,14,1)) + quad_form(p_2,psi_k(:,:,14,1)) + quad_form(p_3,psi_k(:,:,14,1)) + ...
quad_form(p_4,psi_k(:,:,14,2)) + quad_form(p_5,psi_k(:,:,14,2)) + quad_form(p_6,psi_k(:,:,14,2)) + ...
quad_form(p_7,psi_k(:,:,14,3)) + quad_form(p_8,psi_k(:,:,14,3)) + quad_form(p_9,psi_k(:,:,14,3)) + ...
P_n_eff(14)*t_k(14) - 2*real(f_k(:,14)'*p_7) + u_k(14) - v_k(14);

X_15 = quad_form(p_1,psi_k(:,:,15,1)) + quad_form(p_2,psi_k(:,:,15,1)) + quad_form(p_3,psi_k(:,:,15,1)) + ...
quad_form(p_4,psi_k(:,:,15,2)) + quad_form(p_5,psi_k(:,:,15,2)) + quad_form(p_6,psi_k(:,:,15,2)) + ...
quad_form(p_7,psi_k(:,:,15,3)) + quad_form(p_8,psi_k(:,:,15,3)) + quad_form(p_9,psi_k(:,:,15,3)) + ...
P_n_eff(15)*t_k(15) - 2*real(f_k(:,15)'*p_8) + u_k(15) - v_k(15);

X_16 = quad_form(p_1,psi_k(:,:,16,1)) + quad_form(p_2,psi_k(:,:,16,1)) + quad_form(p_3,psi_k(:,:,16,1)) + ...
quad_form(p_4,psi_k(:,:,16,2)) + quad_form(p_5,psi_k(:,:,16,2)) + quad_form(p_6,psi_k(:,:,16,2)) + ...
quad_form(p_7,psi_k(:,:,16,3)) + quad_form(p_8,psi_k(:,:,16,3)) + quad_form(p_9,psi_k(:,:,16,3)) + ...
P_n_eff(16)*t_k(16) - 2*real(f_k(:,16)'*p_8) + u_k(16) - v_k(16);

X_17 = quad_form(p_1,psi_k(:,:,17,1)) + quad_form(p_2,psi_k(:,:,17,1)) + quad_form(p_3,psi_k(:,:,17,1)) + ...
quad_form(p_4,psi_k(:,:,17,2)) + quad_form(p_5,psi_k(:,:,17,2)) + quad_form(p_6,psi_k(:,:,17,2)) + ...
quad_form(p_7,psi_k(:,:,17,3)) + quad_form(p_8,psi_k(:,:,17,3)) + quad_form(p_9,psi_k(:,:,17,3)) + ...
P_n_eff(17)*t_k(17) - 2*real(f_k(:,17)'*p_9) + u_k(17) - v_k(17);

X_18 = quad_form(p_1,psi_k(:,:,18,1)) + quad_form(p_2,psi_k(:,:,18,1)) + quad_form(p_3,psi_k(:,:,18,1)) + ...
quad_form(p_4,psi_k(:,:,18,2)) + quad_form(p_5,psi_k(:,:,18,2)) + quad_form(p_6,psi_k(:,:,18,2)) + ...
quad_form(p_7,psi_k(:,:,18,3)) + quad_form(p_8,psi_k(:,:,18,3)) + quad_form(p_9,psi_k(:,:,18,3)) + ...
P_n_eff(18)*t_k(18) - 2*real(f_k(:,18)'*p_9) + u_k(18) - v_k(18);

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

constraints(28) = PAC - (p_1'*G(:,:,1,1)'*D1*G(:,:,1,1)*p_1 + p_2'*G(:,:,1,1)'*D1*G(:,:,1,1)*p_2 + p_3'*G(:,:,1,1)'*D1*G(:,:,1,1)*p_3 + ... 
p_4'*G(:,:,1,2)'*D1*G(:,:,1,2)*p_4 + p_5'*G(:,:,1,2)'*D1*G(:,:,1,2)*p_5 + p_6'*G(:,:,1,2)'*D1*G(:,:,1,2)*p_6 + ...
p_7'*G(:,:,1,3)'*D1*G(:,:,1,3)*p_7 + p_8'*G(:,:,1,3)'*D1*G(:,:,1,3)*p_8 + p_9'*G(:,:,1,3)'*D1*G(:,:,1,3)*p_9 + 1);

constraints(29) = PAC - (p_1'*G(:,:,1,1)'*D2*G(:,:,1,1)*p_1 + p_2'*G(:,:,1,1)'*D2*G(:,:,1,1)*p_2 + p_3'*G(:,:,1,1)'*D2*G(:,:,1,1)*p_3 + ... 
p_4'*G(:,:,1,2)'*D2*G(:,:,1,2)*p_4 + p_5'*G(:,:,1,2)'*D2*G(:,:,1,2)*p_5 + p_6'*G(:,:,1,2)'*D2*G(:,:,1,2)*p_6 + ...
p_7'*G(:,:,1,3)'*D2*G(:,:,1,3)*p_7 + p_8'*G(:,:,1,3)'*D2*G(:,:,1,3)*p_8 + p_9'*G(:,:,1,3)'*D2*G(:,:,1,3)*p_9 + 1);

constraints(30) = PAC - (p_1'*G(:,:,1,1)'*D3*G(:,:,1,1)*p_1 + p_2'*G(:,:,1,1)'*D3*G(:,:,1,1)*p_2 + p_3'*G(:,:,1,1)'*D3*G(:,:,1,1)*p_3 + ... 
p_4'*G(:,:,1,2)'*D3*G(:,:,1,2)*p_4 + p_5'*G(:,:,1,2)'*D3*G(:,:,1,2)*p_5 + p_6'*G(:,:,1,2)'*D3*G(:,:,1,2)*p_6 + ...
p_7'*G(:,:,1,3)'*D3*G(:,:,1,3)*p_7 + p_8'*G(:,:,1,3)'*D3*G(:,:,1,3)*p_8 + p_9'*G(:,:,1,3)'*D3*G(:,:,1,3)*p_9 + 1);

constraints(31) = PAC - (p_1'*G(:,:,2,1)'*D1*G(:,:,2,1)*p_1 + p_2'*G(:,:,2,1)'*D1*G(:,:,2,1)*p_2 + p_3'*G(:,:,2,1)'*D1*G(:,:,2,1)*p_3 + ... 
p_4'*G(:,:,2,2)'*D1*G(:,:,2,2)*p_4 + p_5'*G(:,:,2,2)'*D1*G(:,:,2,2)*p_5 + p_6'*G(:,:,2,2)'*D1*G(:,:,2,2)*p_6 + ...
p_7'*G(:,:,2,3)'*D1*G(:,:,2,3)*p_7 + p_8'*G(:,:,2,3)'*D1*G(:,:,2,3)*p_8 + p_9'*G(:,:,2,3)'*D1*G(:,:,2,3)*p_9 + 1);

constraints(32) = PAC - (p_1'*G(:,:,2,1)'*D2*G(:,:,2,1)*p_1 + p_2'*G(:,:,2,1)'*D2*G(:,:,2,1)*p_2 + p_3'*G(:,:,2,1)'*D2*G(:,:,2,1)*p_3 + ... 
p_4'*G(:,:,2,2)'*D2*G(:,:,2,2)*p_4 + p_5'*G(:,:,2,2)'*D2*G(:,:,2,2)*p_5 + p_6'*G(:,:,2,2)'*D2*G(:,:,2,2)*p_6 + ...
p_7'*G(:,:,2,3)'*D2*G(:,:,2,3)*p_7 + p_8'*G(:,:,2,3)'*D2*G(:,:,2,3)*p_8 + p_9'*G(:,:,2,3)'*D2*G(:,:,2,3)*p_9 + 1);

constraints(33) = PAC - (p_1'*G(:,:,2,1)'*D3*G(:,:,2,1)*p_1 + p_2'*G(:,:,2,1)'*D3*G(:,:,2,1)*p_2 + p_3'*G(:,:,2,1)'*D3*G(:,:,2,1)*p_3 + ... 
p_4'*G(:,:,2,2)'*D3*G(:,:,2,2)*p_4 + p_5'*G(:,:,2,2)'*D3*G(:,:,2,2)*p_5 + p_6'*G(:,:,2,2)'*D3*G(:,:,2,2)*p_6 + ...
p_7'*G(:,:,2,3)'*D3*G(:,:,2,3)*p_7 + p_8'*G(:,:,2,3)'*D3*G(:,:,2,3)*p_8 + p_9'*G(:,:,2,3)'*D3*G(:,:,2,3)*p_9 + 1);

constraints(34) = PAC - (p_1'*G(:,:,3,1)'*D1*G(:,:,3,1)*p_1 + p_2'*G(:,:,3,1)'*D1*G(:,:,3,1)*p_2 + p_3'*G(:,:,3,1)'*D1*G(:,:,3,1)*p_3 + ... 
p_4'*G(:,:,3,2)'*D1*G(:,:,3,2)*p_4 + p_5'*G(:,:,3,2)'*D1*G(:,:,3,2)*p_5 + p_6'*G(:,:,3,2)'*D1*G(:,:,3,2)*p_6 + ...
p_7'*G(:,:,3,3)'*D1*G(:,:,3,3)*p_7 + p_8'*G(:,:,3,3)'*D1*G(:,:,3,3)*p_8 + p_9'*G(:,:,3,3)'*D1*G(:,:,3,3)*p_9 + 1);

constraints(35) = PAC - (p_1'*G(:,:,3,1)'*D2*G(:,:,3,1)*p_1 + p_2'*G(:,:,3,1)'*D2*G(:,:,3,1)*p_2 + p_3'*G(:,:,3,1)'*D2*G(:,:,3,1)*p_3 + ... 
p_4'*G(:,:,3,2)'*D2*G(:,:,3,2)*p_4 + p_5'*G(:,:,3,2)'*D2*G(:,:,3,2)*p_5 + p_6'*G(:,:,3,2)'*D2*G(:,:,3,2)*p_6 + ...
p_7'*G(:,:,3,3)'*D2*G(:,:,3,3)*p_7 + p_8'*G(:,:,3,3)'*D2*G(:,:,3,3)*p_8 + p_9'*G(:,:,3,3)'*D2*G(:,:,3,3)*p_9 + 1);

constraints(36) = PAC - (p_1'*G(:,:,3,1)'*D3*G(:,:,3,1)*p_1 + p_2'*G(:,:,3,1)'*D3*G(:,:,3,1)*p_2 + p_3'*G(:,:,3,1)'*D3*G(:,:,3,1)*p_3 + ... 
p_4'*G(:,:,3,2)'*D3*G(:,:,3,2)*p_4 + p_5'*G(:,:,3,2)'*D3*G(:,:,3,2)*p_5 + p_6'*G(:,:,3,2)'*D3*G(:,:,3,2)*p_6 + ...
p_7'*G(:,:,3,3)'*D3*G(:,:,3,3)*p_7 + p_8'*G(:,:,3,3)'*D3*G(:,:,3,3)*p_8 + p_9'*G(:,:,3,3)'*D3*G(:,:,3,3)*p_9 + 1);


subject to
    constraints >= zeros(size(constraints))

cvx_end

%Result
rate = object_func;
