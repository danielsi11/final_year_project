function [rate,p_1,p_2,p_3,p_4,p_5,p_6,p_7,p_8,p_9] = NoRS_imperfect_optimisation(H_est,PAC,t_k,psi_k,f_k,v_k,u_k)
    
[N,K] = size(H_est); 

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
X_1 = quad_form(p_1,psi_k(:,:,1)) + quad_form(p_2,psi_k(:,:,1)) + quad_form(p_3,psi_k(:,:,1)) + ...
quad_form(p_4,psi_k(:,:,1)) + quad_form(p_5,psi_k(:,:,1)) + quad_form(p_6,psi_k(:,:,1)) + ...
quad_form(p_7,psi_k(:,:,1)) + quad_form(p_8,psi_k(:,:,1)) + quad_form(p_9,psi_k(:,:,1)) + ...
t_k(1) - 2*real(f_k(:,1)'*p_1) + u_k(1) - v_k(1);

X_2 = quad_form(p_1,psi_k(:,:,2)) + quad_form(p_2,psi_k(:,:,2)) + quad_form(p_3,psi_k(:,:,2)) + ...
quad_form(p_4,psi_k(:,:,2)) + quad_form(p_5,psi_k(:,:,2)) + quad_form(p_6,psi_k(:,:,2)) + ...
quad_form(p_7,psi_k(:,:,2)) + quad_form(p_8,psi_k(:,:,2)) + quad_form(p_9,psi_k(:,:,2)) + ...
t_k(2) - 2*real(f_k(:,2)'*p_1) + u_k(2) - v_k(2);

X_3 = quad_form(p_1,psi_k(:,:,3)) + quad_form(p_2,psi_k(:,:,3)) + quad_form(p_3,psi_k(:,:,3)) + ...
quad_form(p_4,psi_k(:,:,3)) + quad_form(p_5,psi_k(:,:,3)) + quad_form(p_6,psi_k(:,:,3)) + ...
quad_form(p_7,psi_k(:,:,3)) + quad_form(p_8,psi_k(:,:,3)) + quad_form(p_9,psi_k(:,:,3)) + ...
t_k(3) - 2*real(f_k(:,3)'*p_2) + u_k(3) - v_k(3);

X_4 = quad_form(p_1,psi_k(:,:,4)) + quad_form(p_2,psi_k(:,:,4)) + quad_form(p_3,psi_k(:,:,4)) + ...
quad_form(p_4,psi_k(:,:,4)) + quad_form(p_5,psi_k(:,:,4)) + quad_form(p_6,psi_k(:,:,4)) + ...
quad_form(p_7,psi_k(:,:,4)) + quad_form(p_8,psi_k(:,:,4)) + quad_form(p_9,psi_k(:,:,4)) + ...
t_k(4) - 2*real(f_k(:,4)'*p_2) + u_k(4) - v_k(4);

X_5 = quad_form(p_1,psi_k(:,:,5)) + quad_form(p_2,psi_k(:,:,5)) + quad_form(p_3,psi_k(:,:,5)) + ...
quad_form(p_4,psi_k(:,:,5)) + quad_form(p_5,psi_k(:,:,5)) + quad_form(p_6,psi_k(:,:,5)) + ...
quad_form(p_7,psi_k(:,:,5)) + quad_form(p_8,psi_k(:,:,5)) + quad_form(p_9,psi_k(:,:,5)) + ...
t_k(5) - 2*real(f_k(:,5)'*p_3) + u_k(5) - v_k(5);

X_6 = quad_form(p_1,psi_k(:,:,6)) + quad_form(p_2,psi_k(:,:,6)) + quad_form(p_3,psi_k(:,:,6)) + ...
quad_form(p_4,psi_k(:,:,6)) + quad_form(p_5,psi_k(:,:,6)) + quad_form(p_6,psi_k(:,:,6)) + ...
quad_form(p_7,psi_k(:,:,6)) + quad_form(p_8,psi_k(:,:,6)) + quad_form(p_9,psi_k(:,:,6)) + ...
t_k(6) - 2*real(f_k(:,6)'*p_3) + u_k(6) - v_k(6);

X_7 = quad_form(p_1,psi_k(:,:,7)) + quad_form(p_2,psi_k(:,:,7)) + quad_form(p_3,psi_k(:,:,7)) + ...
quad_form(p_4,psi_k(:,:,7)) + quad_form(p_5,psi_k(:,:,7)) + quad_form(p_6,psi_k(:,:,7)) + ...
quad_form(p_7,psi_k(:,:,7)) + quad_form(p_8,psi_k(:,:,7)) + quad_form(p_9,psi_k(:,:,7)) + ...
t_k(7) - 2*real(f_k(:,7)'*p_4) + u_k(7) - v_k(7);

X_8 = quad_form(p_1,psi_k(:,:,8)) + quad_form(p_2,psi_k(:,:,8)) + quad_form(p_3,psi_k(:,:,8)) + ...
quad_form(p_4,psi_k(:,:,8)) + quad_form(p_5,psi_k(:,:,8)) + quad_form(p_6,psi_k(:,:,8)) + ...
quad_form(p_7,psi_k(:,:,8)) + quad_form(p_8,psi_k(:,:,8)) + quad_form(p_9,psi_k(:,:,8)) + ...
t_k(8) - 2*real(f_k(:,8)'*p_4) + u_k(8) - v_k(8);

X_9 = quad_form(p_1,psi_k(:,:,9)) + quad_form(p_2,psi_k(:,:,9)) + quad_form(p_3,psi_k(:,:,9)) + ...
quad_form(p_4,psi_k(:,:,9)) + quad_form(p_5,psi_k(:,:,9)) + quad_form(p_6,psi_k(:,:,9)) + ...
quad_form(p_7,psi_k(:,:,9)) + quad_form(p_8,psi_k(:,:,9)) + quad_form(p_9,psi_k(:,:,9)) + ...
t_k(9) - 2*real(f_k(:,9)'*p_5) + u_k(9) - v_k(9);

X_10 = quad_form(p_1,psi_k(:,:,10)) + quad_form(p_2,psi_k(:,:,10)) + quad_form(p_3,psi_k(:,:,10)) + ...
quad_form(p_4,psi_k(:,:,10)) + quad_form(p_5,psi_k(:,:,10)) + quad_form(p_6,psi_k(:,:,10)) + ...
quad_form(p_7,psi_k(:,:,10)) + quad_form(p_8,psi_k(:,:,10)) + quad_form(p_9,psi_k(:,:,10)) + ...
t_k(10) - 2*real(f_k(:,10)'*p_5) + u_k(10) - v_k(10);

X_11 = quad_form(p_1,psi_k(:,:,11)) + quad_form(p_2,psi_k(:,:,11)) + quad_form(p_3,psi_k(:,:,11)) + ...
quad_form(p_4,psi_k(:,:,11)) + quad_form(p_5,psi_k(:,:,11)) + quad_form(p_6,psi_k(:,:,11)) + ...
quad_form(p_7,psi_k(:,:,11)) + quad_form(p_8,psi_k(:,:,11)) + quad_form(p_9,psi_k(:,:,11)) + ...
t_k(11) - 2*real(f_k(:,11)'*p_6) + u_k(11) - v_k(11);

X_12 = quad_form(p_1,psi_k(:,:,12)) + quad_form(p_2,psi_k(:,:,12)) + quad_form(p_3,psi_k(:,:,12)) + ...
quad_form(p_4,psi_k(:,:,12)) + quad_form(p_5,psi_k(:,:,12)) + quad_form(p_6,psi_k(:,:,12)) + ...
quad_form(p_7,psi_k(:,:,12)) + quad_form(p_8,psi_k(:,:,12)) + quad_form(p_9,psi_k(:,:,12)) + ...
t_k(12) - 2*real(f_k(:,12)'*p_6) + u_k(12) - v_k(12);

X_13 = quad_form(p_1,psi_k(:,:,13)) + quad_form(p_2,psi_k(:,:,13)) + quad_form(p_3,psi_k(:,:,13)) + ...
quad_form(p_4,psi_k(:,:,13)) + quad_form(p_5,psi_k(:,:,13)) + quad_form(p_6,psi_k(:,:,13)) + ...
quad_form(p_7,psi_k(:,:,13)) + quad_form(p_8,psi_k(:,:,13)) + quad_form(p_9,psi_k(:,:,13)) + ...
t_k(13) - 2*real(f_k(:,13)'*p_7) + u_k(13) - v_k(13);

X_14 = quad_form(p_1,psi_k(:,:,14)) + quad_form(p_2,psi_k(:,:,14)) + quad_form(p_3,psi_k(:,:,14)) + ...
quad_form(p_4,psi_k(:,:,14)) + quad_form(p_5,psi_k(:,:,14)) + quad_form(p_6,psi_k(:,:,14)) + ...
quad_form(p_7,psi_k(:,:,14)) + quad_form(p_8,psi_k(:,:,14)) + quad_form(p_9,psi_k(:,:,14)) + ...
t_k(14) - 2*real(f_k(:,14)'*p_7) + u_k(14) - v_k(14);

X_15 = quad_form(p_1,psi_k(:,:,15)) + quad_form(p_2,psi_k(:,:,15)) + quad_form(p_3,psi_k(:,:,15)) + ...
quad_form(p_4,psi_k(:,:,15)) + quad_form(p_5,psi_k(:,:,15)) + quad_form(p_6,psi_k(:,:,15)) + ...
quad_form(p_7,psi_k(:,:,15)) + quad_form(p_8,psi_k(:,:,15)) + quad_form(p_9,psi_k(:,:,15)) + ...
t_k(15) - 2*real(f_k(:,15)'*p_8) + u_k(15) - v_k(15);

X_16 = quad_form(p_1,psi_k(:,:,16)) + quad_form(p_2,psi_k(:,:,16)) + quad_form(p_3,psi_k(:,:,16)) + ...
quad_form(p_4,psi_k(:,:,16)) + quad_form(p_5,psi_k(:,:,16)) + quad_form(p_6,psi_k(:,:,16)) + ...
quad_form(p_7,psi_k(:,:,16)) + quad_form(p_8,psi_k(:,:,16)) + quad_form(p_9,psi_k(:,:,16)) + ...
t_k(16) - 2*real(f_k(:,16)'*p_8) + u_k(16) - v_k(16);

X_17 = quad_form(p_1,psi_k(:,:,17)) + quad_form(p_2,psi_k(:,:,17)) + quad_form(p_3,psi_k(:,:,17)) + ...
quad_form(p_4,psi_k(:,:,17)) + quad_form(p_5,psi_k(:,:,17)) + quad_form(p_6,psi_k(:,:,17)) + ...
quad_form(p_7,psi_k(:,:,17)) + quad_form(p_8,psi_k(:,:,17)) + quad_form(p_9,psi_k(:,:,17)) + ...
t_k(17) - 2*real(f_k(:,17)'*p_9) + u_k(17) - v_k(17);

X_18 = quad_form(p_1,psi_k(:,:,18)) + quad_form(p_2,psi_k(:,:,18)) + quad_form(p_3,psi_k(:,:,18)) + ...
quad_form(p_4,psi_k(:,:,18)) + quad_form(p_5,psi_k(:,:,18)) + quad_form(p_6,psi_k(:,:,18)) + ...
quad_form(p_7,psi_k(:,:,18)) + quad_form(p_8,psi_k(:,:,18)) + quad_form(p_9,psi_k(:,:,18)) + ...
t_k(18) - 2*real(f_k(:,18)'*p_9) + u_k(18) - v_k(18);

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
constraints(28) = PAC - (p_1'*D1*p_1 + p_2'*D1*p_2 + p_3'*D1*p_3 +...
p_4'*D1*p_4 + p_5'*D1*p_5 + p_6'*D1*p_6 + p_7'*D1*p_7 + p_8'*D1*p_8 + p_9'*D1*p_9);
constraints(29) = PAC - (p_1'*D2*p_1 + p_2'*D2*p_2 + p_3'*D2*p_3 +...
p_4'*D2*p_4 + p_5'*D2*p_5 + p_6'*D2*p_6 + p_7'*D2*p_7 + p_8'*D2*p_8 + p_9'*D2*p_9);
constraints(30) = PAC - (p_1'*D3*p_1 + p_2'*D3*p_2 + p_3'*D3*p_3 +...
p_4'*D3*p_4 + p_5'*D3*p_5 + p_6'*D3*p_6 + p_7'*D3*p_7 + p_8'*D3*p_8 + p_9'*D3*p_9);
constraints(31) = PAC - (p_1'*D4*p_1 + p_2'*D4*p_2 + p_3'*D4*p_3 +...
p_4'*D4*p_4 + p_5'*D4*p_5 + p_6'*D4*p_6 + p_7'*D4*p_7 + p_8'*D4*p_8 + p_9'*D4*p_9);
constraints(32) = PAC - (p_1'*D5*p_1 + p_2'*D5*p_2 + p_3'*D5*p_3 +...
p_4'*D5*p_4 + p_5'*D5*p_5 + p_6'*D5*p_6 + p_7'*D5*p_7 + p_8'*D5*p_8 + p_9'*D5*p_9);
constraints(33) = PAC - (p_1'*D6*p_1 + p_2'*D6*p_2 + p_3'*D6*p_3 +...
p_4'*D6*p_4 + p_5'*D6*p_5 + p_6'*D6*p_6 + p_7'*D6*p_7 + p_8'*D6*p_8 + p_9'*D6*p_9);
constraints(34) = PAC - (p_1'*D7*p_1 + p_2'*D7*p_2 + p_3'*D7*p_3 +...
p_4'*D7*p_4 + p_5'*D7*p_5 + p_6'*D7*p_6 + p_7'*D7*p_7 + p_8'*D7*p_8 + p_9'*D7*p_9);
constraints(35) = PAC - (p_1'*D8*p_1 + p_2'*D8*p_2 + p_3'*D8*p_3 +...
p_4'*D8*p_4 + p_5'*D8*p_5 + p_6'*D8*p_6 + p_7'*D8*p_7 + p_8'*D8*p_8 + p_9'*D8*p_9);
constraints(36) = PAC - (p_1'*D9*p_1 + p_2'*D9*p_2 + p_3'*D9*p_3 +...
p_4'*D9*p_4 + p_5'*D9*p_5 + p_6'*D9*p_6 + p_7'*D9*p_7 + p_8'*D9*p_8 + p_9'*D9*p_9);

subject to
    constraints >= zeros(size(constraints))

cvx_end

%Result
rate = object_func;