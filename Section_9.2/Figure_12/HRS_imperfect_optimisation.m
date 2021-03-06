function [rate,p_sc,p_c1,p_c2,p_c3,p_1,p_2,p_3,p_4,p_5,p_6,p_7,p_8,p_9] = ...
HRS_imperfect_optimisation(H_est,PAC,t_k,t_ck,t_sck,psi_k,psi_ck,psi_sck,f_k,f_ck,f_sck,v_k,v_ck,v_sck,u_k,u_ck,u_sck)
    
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
variable p_c1(N,1) complex
variable p_c2(N,1) complex
variable p_c3(N,1) complex
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

expression constraints(1,90);

%Rate-WMMSE Relationship
%Private Stream
X_1 = quad_form(p_c2,psi_k(:,:,1)) + quad_form(p_c3,psi_k(:,:,1)) + ...
quad_form(p_1,psi_k(:,:,1)) + quad_form(p_2,psi_k(:,:,1)) + quad_form(p_3,psi_k(:,:,1)) + ...
quad_form(p_4,psi_k(:,:,1)) + quad_form(p_5,psi_k(:,:,1)) + quad_form(p_6,psi_k(:,:,1)) + ...
quad_form(p_7,psi_k(:,:,1)) + quad_form(p_8,psi_k(:,:,1)) + quad_form(p_9,psi_k(:,:,1)) + ...
t_k(1) - 2*real(f_k(:,1)'*p_1) + u_k(1) - v_k(1);

X_2 = quad_form(p_c2,psi_k(:,:,2)) + quad_form(p_c3,psi_k(:,:,2)) + ...
quad_form(p_1,psi_k(:,:,2)) + quad_form(p_2,psi_k(:,:,2)) + quad_form(p_3,psi_k(:,:,2)) + ...
quad_form(p_4,psi_k(:,:,2)) + quad_form(p_5,psi_k(:,:,2)) + quad_form(p_6,psi_k(:,:,2)) + ...
quad_form(p_7,psi_k(:,:,2)) + quad_form(p_8,psi_k(:,:,2)) + quad_form(p_9,psi_k(:,:,2)) + ...
t_k(2) - 2*real(f_k(:,2)'*p_1) + u_k(2) - v_k(2);

X_3 = quad_form(p_c2,psi_k(:,:,3)) + quad_form(p_c3,psi_k(:,:,3)) + ...
quad_form(p_1,psi_k(:,:,3)) + quad_form(p_2,psi_k(:,:,3)) + quad_form(p_3,psi_k(:,:,3)) + ...
quad_form(p_4,psi_k(:,:,3)) + quad_form(p_5,psi_k(:,:,3)) + quad_form(p_6,psi_k(:,:,3)) + ...
quad_form(p_7,psi_k(:,:,3)) + quad_form(p_8,psi_k(:,:,3)) + quad_form(p_9,psi_k(:,:,3)) + ...
t_k(3) - 2*real(f_k(:,3)'*p_2) + u_k(3) - v_k(3);

X_4 = quad_form(p_c2,psi_k(:,:,4)) + quad_form(p_c3,psi_k(:,:,4)) + ...
quad_form(p_1,psi_k(:,:,4)) + quad_form(p_2,psi_k(:,:,4)) + quad_form(p_3,psi_k(:,:,4)) + ...
quad_form(p_4,psi_k(:,:,4)) + quad_form(p_5,psi_k(:,:,4)) + quad_form(p_6,psi_k(:,:,4)) + ...
quad_form(p_7,psi_k(:,:,4)) + quad_form(p_8,psi_k(:,:,4)) + quad_form(p_9,psi_k(:,:,4)) + ...
t_k(4) - 2*real(f_k(:,4)'*p_2) + u_k(4) - v_k(4);

X_5 = quad_form(p_c2,psi_k(:,:,5)) + quad_form(p_c3,psi_k(:,:,5)) + ...
quad_form(p_1,psi_k(:,:,5)) + quad_form(p_2,psi_k(:,:,5)) + quad_form(p_3,psi_k(:,:,5)) + ...
quad_form(p_4,psi_k(:,:,5)) + quad_form(p_5,psi_k(:,:,5)) + quad_form(p_6,psi_k(:,:,5)) + ...
quad_form(p_7,psi_k(:,:,5)) + quad_form(p_8,psi_k(:,:,5)) + quad_form(p_9,psi_k(:,:,5)) + ...
t_k(5) - 2*real(f_k(:,5)'*p_3) + u_k(5) - v_k(5);

X_6 = quad_form(p_c2,psi_k(:,:,6)) + quad_form(p_c3,psi_k(:,:,6)) + ...
quad_form(p_1,psi_k(:,:,6)) + quad_form(p_2,psi_k(:,:,6)) + quad_form(p_3,psi_k(:,:,6)) + ...
quad_form(p_4,psi_k(:,:,6)) + quad_form(p_5,psi_k(:,:,6)) + quad_form(p_6,psi_k(:,:,6)) + ...
quad_form(p_7,psi_k(:,:,6)) + quad_form(p_8,psi_k(:,:,6)) + quad_form(p_9,psi_k(:,:,6)) + ...
t_k(6) - 2*real(f_k(:,6)'*p_3) + u_k(6) - v_k(6);

X_7 = quad_form(p_c1,psi_k(:,:,7)) + quad_form(p_c3,psi_k(:,:,7)) + ...
quad_form(p_1,psi_k(:,:,7)) + quad_form(p_2,psi_k(:,:,7)) + quad_form(p_3,psi_k(:,:,7)) + ...
quad_form(p_4,psi_k(:,:,7)) + quad_form(p_5,psi_k(:,:,7)) + quad_form(p_6,psi_k(:,:,7)) + ...
quad_form(p_7,psi_k(:,:,7)) + quad_form(p_8,psi_k(:,:,7)) + quad_form(p_9,psi_k(:,:,7)) + ...
t_k(7) - 2*real(f_k(:,7)'*p_4) + u_k(7) - v_k(7);

X_8 = quad_form(p_c1,psi_k(:,:,8)) + quad_form(p_c3,psi_k(:,:,8)) + ...
quad_form(p_1,psi_k(:,:,8)) + quad_form(p_2,psi_k(:,:,8)) + quad_form(p_3,psi_k(:,:,8)) + ...
quad_form(p_4,psi_k(:,:,8)) + quad_form(p_5,psi_k(:,:,8)) + quad_form(p_6,psi_k(:,:,8)) + ...
quad_form(p_7,psi_k(:,:,8)) + quad_form(p_8,psi_k(:,:,8)) + quad_form(p_9,psi_k(:,:,8)) + ...
t_k(8) - 2*real(f_k(:,8)'*p_4) + u_k(8) - v_k(8);

X_9 = quad_form(p_c1,psi_k(:,:,9)) + quad_form(p_c3,psi_k(:,:,9)) + ...
quad_form(p_1,psi_k(:,:,9)) + quad_form(p_2,psi_k(:,:,9)) + quad_form(p_3,psi_k(:,:,9)) + ...
quad_form(p_4,psi_k(:,:,9)) + quad_form(p_5,psi_k(:,:,9)) + quad_form(p_6,psi_k(:,:,9)) + ...
quad_form(p_7,psi_k(:,:,9)) + quad_form(p_8,psi_k(:,:,9)) + quad_form(p_9,psi_k(:,:,9)) + ...
t_k(9) - 2*real(f_k(:,9)'*p_5) + u_k(9) - v_k(9);

X_10 = quad_form(p_c1,psi_k(:,:,10)) + quad_form(p_c3,psi_k(:,:,10)) + ...
quad_form(p_1,psi_k(:,:,10)) + quad_form(p_2,psi_k(:,:,10)) + quad_form(p_3,psi_k(:,:,10)) + ...
quad_form(p_4,psi_k(:,:,10)) + quad_form(p_5,psi_k(:,:,10)) + quad_form(p_6,psi_k(:,:,10)) + ...
quad_form(p_7,psi_k(:,:,10)) + quad_form(p_8,psi_k(:,:,10)) + quad_form(p_9,psi_k(:,:,10)) + ...
t_k(10) - 2*real(f_k(:,10)'*p_5) + u_k(10) - v_k(10);

X_11 = quad_form(p_c1,psi_k(:,:,11)) + quad_form(p_c3,psi_k(:,:,11)) + ...
quad_form(p_1,psi_k(:,:,11)) + quad_form(p_2,psi_k(:,:,11)) + quad_form(p_3,psi_k(:,:,11)) + ...
quad_form(p_4,psi_k(:,:,11)) + quad_form(p_5,psi_k(:,:,11)) + quad_form(p_6,psi_k(:,:,11)) + ...
quad_form(p_7,psi_k(:,:,11)) + quad_form(p_8,psi_k(:,:,11)) + quad_form(p_9,psi_k(:,:,11)) + ...
t_k(11) - 2*real(f_k(:,11)'*p_6) + u_k(11) - v_k(11);

X_12 = quad_form(p_c1,psi_k(:,:,12)) + quad_form(p_c3,psi_k(:,:,12)) + ...
quad_form(p_1,psi_k(:,:,12)) + quad_form(p_2,psi_k(:,:,12)) + quad_form(p_3,psi_k(:,:,12)) + ...
quad_form(p_4,psi_k(:,:,12)) + quad_form(p_5,psi_k(:,:,12)) + quad_form(p_6,psi_k(:,:,12)) + ...
quad_form(p_7,psi_k(:,:,12)) + quad_form(p_8,psi_k(:,:,12)) + quad_form(p_9,psi_k(:,:,12)) + ...
t_k(12) - 2*real(f_k(:,12)'*p_6) + u_k(12) - v_k(12);

X_13 = quad_form(p_c1,psi_k(:,:,13)) + quad_form(p_c2,psi_k(:,:,13)) + ...
quad_form(p_1,psi_k(:,:,13)) + quad_form(p_2,psi_k(:,:,13)) + quad_form(p_3,psi_k(:,:,13)) + ...
quad_form(p_4,psi_k(:,:,13)) + quad_form(p_5,psi_k(:,:,13)) + quad_form(p_6,psi_k(:,:,13)) + ...
quad_form(p_7,psi_k(:,:,13)) + quad_form(p_8,psi_k(:,:,13)) + quad_form(p_9,psi_k(:,:,13)) + ...
t_k(13) - 2*real(f_k(:,13)'*p_7) + u_k(13) - v_k(13);

X_14 = quad_form(p_c1,psi_k(:,:,14)) + quad_form(p_c2,psi_k(:,:,14)) + ...
quad_form(p_1,psi_k(:,:,14)) + quad_form(p_2,psi_k(:,:,14)) + quad_form(p_3,psi_k(:,:,14)) + ...
quad_form(p_4,psi_k(:,:,14)) + quad_form(p_5,psi_k(:,:,14)) + quad_form(p_6,psi_k(:,:,14)) + ...
quad_form(p_7,psi_k(:,:,14)) + quad_form(p_8,psi_k(:,:,14)) + quad_form(p_9,psi_k(:,:,14)) + ...
t_k(14) - 2*real(f_k(:,14)'*p_7) + u_k(14) - v_k(14);

X_15 = quad_form(p_c1,psi_k(:,:,15)) + quad_form(p_c2,psi_k(:,:,15)) + ...
quad_form(p_1,psi_k(:,:,15)) + quad_form(p_2,psi_k(:,:,15)) + quad_form(p_3,psi_k(:,:,15)) + ...
quad_form(p_4,psi_k(:,:,15)) + quad_form(p_5,psi_k(:,:,15)) + quad_form(p_6,psi_k(:,:,15)) + ...
quad_form(p_7,psi_k(:,:,15)) + quad_form(p_8,psi_k(:,:,15)) + quad_form(p_9,psi_k(:,:,15)) + ...
t_k(15) - 2*real(f_k(:,15)'*p_8) + u_k(15) - v_k(15);

X_16 = quad_form(p_c1,psi_k(:,:,16)) + quad_form(p_c2,psi_k(:,:,16)) + ...
quad_form(p_1,psi_k(:,:,16)) + quad_form(p_2,psi_k(:,:,16)) + quad_form(p_3,psi_k(:,:,16)) + ...
quad_form(p_4,psi_k(:,:,16)) + quad_form(p_5,psi_k(:,:,16)) + quad_form(p_6,psi_k(:,:,16)) + ...
quad_form(p_7,psi_k(:,:,16)) + quad_form(p_8,psi_k(:,:,16)) + quad_form(p_9,psi_k(:,:,16)) + ...
t_k(16) - 2*real(f_k(:,16)'*p_8) + u_k(16) - v_k(16);

X_17 = quad_form(p_c1,psi_k(:,:,17)) + quad_form(p_c2,psi_k(:,:,17)) + ...
quad_form(p_1,psi_k(:,:,17)) + quad_form(p_2,psi_k(:,:,17)) + quad_form(p_3,psi_k(:,:,17)) + ...
quad_form(p_4,psi_k(:,:,17)) + quad_form(p_5,psi_k(:,:,17)) + quad_form(p_6,psi_k(:,:,17)) + ...
quad_form(p_7,psi_k(:,:,17)) + quad_form(p_8,psi_k(:,:,17)) + quad_form(p_9,psi_k(:,:,17)) + ...
t_k(17) - 2*real(f_k(:,17)'*p_9) + u_k(17) - v_k(17);

X_18 = quad_form(p_c1,psi_k(:,:,18)) + quad_form(p_c2,psi_k(:,:,18)) + ...
quad_form(p_1,psi_k(:,:,18)) + quad_form(p_2,psi_k(:,:,18)) + quad_form(p_3,psi_k(:,:,18)) + ...
quad_form(p_4,psi_k(:,:,18)) + quad_form(p_5,psi_k(:,:,18)) + quad_form(p_6,psi_k(:,:,18)) + ...
quad_form(p_7,psi_k(:,:,18)) + quad_form(p_8,psi_k(:,:,18)) + quad_form(p_9,psi_k(:,:,18)) + ...
t_k(18) - 2*real(f_k(:,18)'*p_9) + u_k(18) - v_k(18);

%Local Common Stream
X_c1 = quad_form(p_c1,psi_ck(:,:,1)) + quad_form(p_c2,psi_ck(:,:,1)) + ...
quad_form(p_c3,psi_ck(:,:,1)) + quad_form(p_1,psi_ck(:,:,1)) + quad_form(p_2,psi_ck(:,:,1)) + ...
quad_form(p_3,psi_ck(:,:,1)) + quad_form(p_4,psi_ck(:,:,1)) + quad_form(p_5,psi_ck(:,:,1)) + ...
quad_form(p_6,psi_ck(:,:,1)) + quad_form(p_7,psi_ck(:,:,1)) + quad_form(p_8,psi_ck(:,:,1)) + ...
quad_form(p_9,psi_ck(:,:,1)) + t_ck(1) - 2*real(f_ck(:,1)'*p_c1) + u_ck(1) - v_ck(1);

X_c2 = quad_form(p_c1,psi_ck(:,:,2)) + quad_form(p_c2,psi_ck(:,:,2)) + ...
quad_form(p_c3,psi_ck(:,:,2)) + quad_form(p_1,psi_ck(:,:,2)) + quad_form(p_2,psi_ck(:,:,2)) + ...
quad_form(p_3,psi_ck(:,:,2)) + quad_form(p_4,psi_ck(:,:,2)) + quad_form(p_5,psi_ck(:,:,2)) + ...
quad_form(p_6,psi_ck(:,:,2)) + quad_form(p_7,psi_ck(:,:,2)) + quad_form(p_8,psi_ck(:,:,2)) + ...
quad_form(p_9,psi_ck(:,:,2)) + t_ck(2) - 2*real(f_ck(:,2)'*p_c1) + u_ck(2) - v_ck(2);

X_c3 = quad_form(p_c1,psi_ck(:,:,3)) + quad_form(p_c2,psi_ck(:,:,3)) + ...
quad_form(p_c3,psi_ck(:,:,3)) + quad_form(p_1,psi_ck(:,:,3)) + quad_form(p_2,psi_ck(:,:,3)) + ...
quad_form(p_3,psi_ck(:,:,3)) + quad_form(p_4,psi_ck(:,:,3)) + quad_form(p_5,psi_ck(:,:,3)) + ...
quad_form(p_6,psi_ck(:,:,3)) + quad_form(p_7,psi_ck(:,:,3)) + quad_form(p_8,psi_ck(:,:,3)) + ...
quad_form(p_9,psi_ck(:,:,3)) + t_ck(3) - 2*real(f_ck(:,3)'*p_c1) + u_ck(3) - v_ck(3);

X_c4 = quad_form(p_c1,psi_ck(:,:,4)) + quad_form(p_c2,psi_ck(:,:,4)) + ...
quad_form(p_c3,psi_ck(:,:,4)) + quad_form(p_1,psi_ck(:,:,4)) + quad_form(p_2,psi_ck(:,:,4)) + ...
quad_form(p_3,psi_ck(:,:,4)) + quad_form(p_4,psi_ck(:,:,4)) + quad_form(p_5,psi_ck(:,:,4)) + ...
quad_form(p_6,psi_ck(:,:,4)) + quad_form(p_7,psi_ck(:,:,4)) + quad_form(p_8,psi_ck(:,:,4)) + ...
quad_form(p_9,psi_ck(:,:,4)) + t_ck(4) - 2*real(f_ck(:,4)'*p_c1) + u_ck(4) - v_ck(4);

X_c5 = quad_form(p_c1,psi_ck(:,:,5)) + quad_form(p_c2,psi_ck(:,:,5)) + ...
quad_form(p_c3,psi_ck(:,:,5)) + quad_form(p_1,psi_ck(:,:,5)) + quad_form(p_2,psi_ck(:,:,5)) + ...
quad_form(p_3,psi_ck(:,:,5)) + quad_form(p_4,psi_ck(:,:,5)) + quad_form(p_5,psi_ck(:,:,5)) + ...
quad_form(p_6,psi_ck(:,:,5)) + quad_form(p_7,psi_ck(:,:,5)) + quad_form(p_8,psi_ck(:,:,5)) + ...
quad_form(p_9,psi_ck(:,:,5)) + t_ck(5) - 2*real(f_ck(:,5)'*p_c1) + u_ck(5) - v_ck(5);

X_c6 = quad_form(p_c1,psi_ck(:,:,6)) + quad_form(p_c2,psi_ck(:,:,6)) + ...
quad_form(p_c3,psi_ck(:,:,6)) + quad_form(p_1,psi_ck(:,:,6)) + quad_form(p_2,psi_ck(:,:,6)) + ...
quad_form(p_3,psi_ck(:,:,6)) + quad_form(p_4,psi_ck(:,:,6)) + quad_form(p_5,psi_ck(:,:,6)) + ...
quad_form(p_6,psi_ck(:,:,6)) + quad_form(p_7,psi_ck(:,:,6)) + quad_form(p_8,psi_ck(:,:,6)) + ...
quad_form(p_9,psi_ck(:,:,6)) + t_ck(6) - 2*real(f_ck(:,6)'*p_c1) + u_ck(6) - v_ck(6);

X_c7 = quad_form(p_c1,psi_ck(:,:,7)) + quad_form(p_c2,psi_ck(:,:,7)) + ...
quad_form(p_c3,psi_ck(:,:,7)) + quad_form(p_1,psi_ck(:,:,7)) + quad_form(p_2,psi_ck(:,:,7)) + ...
quad_form(p_3,psi_ck(:,:,7)) + quad_form(p_4,psi_ck(:,:,7)) + quad_form(p_5,psi_ck(:,:,7)) + ...
quad_form(p_6,psi_ck(:,:,7)) + quad_form(p_7,psi_ck(:,:,7)) + quad_form(p_8,psi_ck(:,:,7)) + ...
quad_form(p_9,psi_ck(:,:,7)) + t_ck(7) - 2*real(f_ck(:,7)'*p_c2) + u_ck(7) - v_ck(7);

X_c8 = quad_form(p_c1,psi_ck(:,:,8)) + quad_form(p_c2,psi_ck(:,:,8)) + ...
quad_form(p_c3,psi_ck(:,:,8)) + quad_form(p_1,psi_ck(:,:,8)) + quad_form(p_2,psi_ck(:,:,8)) + ...
quad_form(p_3,psi_ck(:,:,8)) + quad_form(p_4,psi_ck(:,:,8)) + quad_form(p_5,psi_ck(:,:,8)) + ...
quad_form(p_6,psi_ck(:,:,8)) + quad_form(p_7,psi_ck(:,:,8)) + quad_form(p_8,psi_ck(:,:,8)) + ...
quad_form(p_9,psi_ck(:,:,8)) + t_ck(8) - 2*real(f_ck(:,8)'*p_c2) + u_ck(8) - v_ck(8);

X_c9 = quad_form(p_c1,psi_ck(:,:,9)) + quad_form(p_c2,psi_ck(:,:,9)) + ...
quad_form(p_c3,psi_ck(:,:,9)) + quad_form(p_1,psi_ck(:,:,9)) + quad_form(p_2,psi_ck(:,:,9)) + ...
quad_form(p_3,psi_ck(:,:,9)) + quad_form(p_4,psi_ck(:,:,9)) + quad_form(p_5,psi_ck(:,:,9)) + ...
quad_form(p_6,psi_ck(:,:,9)) + quad_form(p_7,psi_ck(:,:,9)) + quad_form(p_8,psi_ck(:,:,9)) + ...
quad_form(p_9,psi_ck(:,:,9)) + t_ck(9) - 2*real(f_ck(:,9)'*p_c2) + u_ck(9) - v_ck(9);

X_c10 = quad_form(p_c1,psi_ck(:,:,10)) + quad_form(p_c2,psi_ck(:,:,10)) + ...
quad_form(p_c3,psi_ck(:,:,10)) + quad_form(p_1,psi_ck(:,:,10)) + quad_form(p_2,psi_ck(:,:,10)) + ...
quad_form(p_3,psi_ck(:,:,10)) + quad_form(p_4,psi_ck(:,:,10)) + quad_form(p_5,psi_ck(:,:,10)) + ...
quad_form(p_6,psi_ck(:,:,10)) + quad_form(p_7,psi_ck(:,:,10)) + quad_form(p_8,psi_ck(:,:,10)) + ...
quad_form(p_9,psi_ck(:,:,10)) + t_ck(10) - 2*real(f_ck(:,10)'*p_c2) + u_ck(10) - v_ck(10);

X_c11 = quad_form(p_c1,psi_ck(:,:,11)) + quad_form(p_c2,psi_ck(:,:,11)) + ...
quad_form(p_c3,psi_ck(:,:,11)) + quad_form(p_1,psi_ck(:,:,11)) + quad_form(p_2,psi_ck(:,:,11)) + ...
quad_form(p_3,psi_ck(:,:,11)) + quad_form(p_4,psi_ck(:,:,11)) + quad_form(p_5,psi_ck(:,:,11)) + ...
quad_form(p_6,psi_ck(:,:,11)) + quad_form(p_7,psi_ck(:,:,11)) + quad_form(p_8,psi_ck(:,:,11)) + ...
quad_form(p_9,psi_ck(:,:,11)) + t_ck(11) - 2*real(f_ck(:,11)'*p_c2) + u_ck(11) - v_ck(11);

X_c12 = quad_form(p_c1,psi_ck(:,:,12)) + quad_form(p_c2,psi_ck(:,:,12)) + ...
quad_form(p_c3,psi_ck(:,:,12)) + quad_form(p_1,psi_ck(:,:,12)) + quad_form(p_2,psi_ck(:,:,12)) + ...
quad_form(p_3,psi_ck(:,:,12)) + quad_form(p_4,psi_ck(:,:,12)) + quad_form(p_5,psi_ck(:,:,12)) + ...
quad_form(p_6,psi_ck(:,:,12)) + quad_form(p_7,psi_ck(:,:,12)) + quad_form(p_8,psi_ck(:,:,12)) + ...
quad_form(p_9,psi_ck(:,:,12)) + t_ck(12) - 2*real(f_ck(:,12)'*p_c2) + u_ck(12) - v_ck(12);

X_c13 = quad_form(p_c1,psi_ck(:,:,13)) + quad_form(p_c2,psi_ck(:,:,13)) + ...
quad_form(p_c3,psi_ck(:,:,13)) + quad_form(p_1,psi_ck(:,:,13)) + quad_form(p_2,psi_ck(:,:,13)) + ...
quad_form(p_3,psi_ck(:,:,13)) + quad_form(p_4,psi_ck(:,:,13)) + quad_form(p_5,psi_ck(:,:,13)) + ...
quad_form(p_6,psi_ck(:,:,13)) + quad_form(p_7,psi_ck(:,:,13)) + quad_form(p_8,psi_ck(:,:,13)) + ...
quad_form(p_9,psi_ck(:,:,13)) + t_ck(13) - 2*real(f_ck(:,13)'*p_c3) + u_ck(13) - v_ck(13);

X_c14 = quad_form(p_c1,psi_ck(:,:,14)) + quad_form(p_c2,psi_ck(:,:,14)) + ...
quad_form(p_c3,psi_ck(:,:,14)) + quad_form(p_1,psi_ck(:,:,14)) + quad_form(p_2,psi_ck(:,:,14)) + ...
quad_form(p_3,psi_ck(:,:,14)) + quad_form(p_4,psi_ck(:,:,14)) + quad_form(p_5,psi_ck(:,:,14)) + ...
quad_form(p_6,psi_ck(:,:,14)) + quad_form(p_7,psi_ck(:,:,14)) + quad_form(p_8,psi_ck(:,:,14)) + ...
quad_form(p_9,psi_ck(:,:,14)) + t_ck(14) - 2*real(f_ck(:,14)'*p_c3) + u_ck(14) - v_ck(14);

X_c15 = quad_form(p_c1,psi_ck(:,:,15)) + quad_form(p_c2,psi_ck(:,:,15)) + ...
quad_form(p_c3,psi_ck(:,:,15)) + quad_form(p_1,psi_ck(:,:,15)) + quad_form(p_2,psi_ck(:,:,15)) + ...
quad_form(p_3,psi_ck(:,:,15)) + quad_form(p_4,psi_ck(:,:,15)) + quad_form(p_5,psi_ck(:,:,15)) + ...
quad_form(p_6,psi_ck(:,:,15)) + quad_form(p_7,psi_ck(:,:,15)) + quad_form(p_8,psi_ck(:,:,15)) + ...
quad_form(p_9,psi_ck(:,:,15)) + t_ck(15) - 2*real(f_ck(:,15)'*p_c3) + u_ck(15) - v_ck(15);

X_c16 = quad_form(p_c1,psi_ck(:,:,16)) + quad_form(p_c2,psi_ck(:,:,16)) + ...
quad_form(p_c3,psi_ck(:,:,16)) + quad_form(p_1,psi_ck(:,:,16)) + quad_form(p_2,psi_ck(:,:,16)) + ...
quad_form(p_3,psi_ck(:,:,16)) + quad_form(p_4,psi_ck(:,:,16)) + quad_form(p_5,psi_ck(:,:,16)) + ...
quad_form(p_6,psi_ck(:,:,16)) + quad_form(p_7,psi_ck(:,:,16)) + quad_form(p_8,psi_ck(:,:,16)) + ...
quad_form(p_9,psi_ck(:,:,16)) + t_ck(16) - 2*real(f_ck(:,16)'*p_c3) + u_ck(16) - v_ck(16);

X_c17 = quad_form(p_c1,psi_ck(:,:,17)) + quad_form(p_c2,psi_ck(:,:,17)) + ...
quad_form(p_c3,psi_ck(:,:,17)) + quad_form(p_1,psi_ck(:,:,17)) + quad_form(p_2,psi_ck(:,:,17)) + ...
quad_form(p_3,psi_ck(:,:,17)) + quad_form(p_4,psi_ck(:,:,17)) + quad_form(p_5,psi_ck(:,:,17)) + ...
quad_form(p_6,psi_ck(:,:,17)) + quad_form(p_7,psi_ck(:,:,17)) + quad_form(p_8,psi_ck(:,:,17)) + ...
quad_form(p_9,psi_ck(:,:,17)) + t_ck(17) - 2*real(f_ck(:,17)'*p_c3) + u_ck(17) - v_ck(17);

X_c18 = quad_form(p_c1,psi_ck(:,:,18)) + quad_form(p_c2,psi_ck(:,:,18)) + ...
quad_form(p_c3,psi_ck(:,:,18)) + quad_form(p_1,psi_ck(:,:,18)) + quad_form(p_2,psi_ck(:,:,18)) + ...
quad_form(p_3,psi_ck(:,:,18)) + quad_form(p_4,psi_ck(:,:,18)) + quad_form(p_5,psi_ck(:,:,18)) + ...
quad_form(p_6,psi_ck(:,:,18)) + quad_form(p_7,psi_ck(:,:,18)) + quad_form(p_8,psi_ck(:,:,18)) + ...
quad_form(p_9,psi_ck(:,:,18)) + t_ck(18) - 2*real(f_ck(:,18)'*p_c3) + u_ck(18) - v_ck(18);

%Global Common Stream
X_sc1 = quad_form(p_sc,psi_sck(:,:,1)) + quad_form(p_c1,psi_sck(:,:,1)) + quad_form(p_c2,psi_sck(:,:,1)) + ...
quad_form(p_c3,psi_sck(:,:,1)) + quad_form(p_1,psi_sck(:,:,1)) + quad_form(p_2,psi_sck(:,:,1)) + ...
quad_form(p_3,psi_sck(:,:,1)) + quad_form(p_4,psi_sck(:,:,1)) + quad_form(p_5,psi_sck(:,:,1)) + ...
quad_form(p_6,psi_sck(:,:,1)) + quad_form(p_7,psi_sck(:,:,1)) + quad_form(p_8,psi_sck(:,:,1)) + ...
quad_form(p_9,psi_sck(:,:,1)) + t_sck(1) - 2*real(f_sck(:,1)'*p_sc) + u_sck(1) - v_sck(1);

X_sc2 = quad_form(p_sc,psi_sck(:,:,2)) + quad_form(p_c1,psi_sck(:,:,2)) + quad_form(p_c2,psi_sck(:,:,2)) + ...
quad_form(p_c3,psi_sck(:,:,2)) + quad_form(p_1,psi_sck(:,:,2)) + quad_form(p_2,psi_sck(:,:,2)) + ...
quad_form(p_3,psi_sck(:,:,2)) + quad_form(p_4,psi_sck(:,:,2)) + quad_form(p_5,psi_sck(:,:,2)) + ...
quad_form(p_6,psi_sck(:,:,2)) + quad_form(p_7,psi_sck(:,:,2)) + quad_form(p_8,psi_sck(:,:,2)) + ...
quad_form(p_9,psi_sck(:,:,2)) + t_sck(2) - 2*real(f_sck(:,2)'*p_sc) + u_sck(2) - v_sck(2);

X_sc3 = quad_form(p_sc,psi_sck(:,:,3)) + quad_form(p_c1,psi_sck(:,:,3)) + quad_form(p_c2,psi_sck(:,:,3)) + ...
quad_form(p_c3,psi_sck(:,:,3)) + quad_form(p_1,psi_sck(:,:,3)) + quad_form(p_2,psi_sck(:,:,3)) + ...
quad_form(p_3,psi_sck(:,:,3)) + quad_form(p_4,psi_sck(:,:,3)) + quad_form(p_5,psi_sck(:,:,3)) + ...
quad_form(p_6,psi_sck(:,:,3)) + quad_form(p_7,psi_sck(:,:,3)) + quad_form(p_8,psi_sck(:,:,3)) + ...
quad_form(p_9,psi_sck(:,:,3)) + t_sck(3) - 2*real(f_sck(:,3)'*p_sc) + u_sck(3) - v_sck(3);

X_sc4 = quad_form(p_sc,psi_sck(:,:,4)) + quad_form(p_c1,psi_sck(:,:,4)) + quad_form(p_c2,psi_sck(:,:,4)) + ...
quad_form(p_c3,psi_sck(:,:,4)) + quad_form(p_1,psi_sck(:,:,4)) + quad_form(p_2,psi_sck(:,:,4)) + ...
quad_form(p_3,psi_sck(:,:,4)) + quad_form(p_4,psi_sck(:,:,4)) + quad_form(p_5,psi_sck(:,:,4)) + ...
quad_form(p_6,psi_sck(:,:,4)) + quad_form(p_7,psi_sck(:,:,4)) + quad_form(p_8,psi_sck(:,:,4)) + ...
quad_form(p_9,psi_sck(:,:,4)) + t_sck(4) - 2*real(f_sck(:,4)'*p_sc) + u_sck(4) - v_sck(4);

X_sc5 = quad_form(p_sc,psi_sck(:,:,5)) + quad_form(p_c1,psi_sck(:,:,5)) + quad_form(p_c2,psi_sck(:,:,5)) + ...
quad_form(p_c3,psi_sck(:,:,5)) + quad_form(p_1,psi_sck(:,:,5)) + quad_form(p_2,psi_sck(:,:,5)) + ...
quad_form(p_3,psi_sck(:,:,5)) + quad_form(p_4,psi_sck(:,:,5)) + quad_form(p_5,psi_sck(:,:,5)) + ...
quad_form(p_6,psi_sck(:,:,5)) + quad_form(p_7,psi_sck(:,:,5)) + quad_form(p_8,psi_sck(:,:,5)) + ...
quad_form(p_9,psi_sck(:,:,5)) + t_sck(5) - 2*real(f_sck(:,5)'*p_sc) + u_sck(5) - v_sck(5);

X_sc6 = quad_form(p_sc,psi_sck(:,:,6)) + quad_form(p_c1,psi_sck(:,:,6)) + quad_form(p_c2,psi_sck(:,:,6)) + ...
quad_form(p_c3,psi_sck(:,:,6)) + quad_form(p_1,psi_sck(:,:,6)) + quad_form(p_2,psi_sck(:,:,6)) + ...
quad_form(p_3,psi_sck(:,:,6)) + quad_form(p_4,psi_sck(:,:,6)) + quad_form(p_5,psi_sck(:,:,6)) + ...
quad_form(p_6,psi_sck(:,:,6)) + quad_form(p_7,psi_sck(:,:,6)) + quad_form(p_8,psi_sck(:,:,6)) + ...
quad_form(p_9,psi_sck(:,:,6)) + t_sck(6) - 2*real(f_sck(:,6)'*p_sc) + u_sck(6) - v_sck(6);

X_sc7 = quad_form(p_sc,psi_sck(:,:,7)) + quad_form(p_c1,psi_sck(:,:,7)) + quad_form(p_c2,psi_sck(:,:,7)) + ...
quad_form(p_c3,psi_sck(:,:,7)) + quad_form(p_1,psi_sck(:,:,7)) + quad_form(p_2,psi_sck(:,:,7)) + ...
quad_form(p_3,psi_sck(:,:,7)) + quad_form(p_4,psi_sck(:,:,7)) + quad_form(p_5,psi_sck(:,:,7)) + ...
quad_form(p_6,psi_sck(:,:,7)) + quad_form(p_7,psi_sck(:,:,7)) + quad_form(p_8,psi_sck(:,:,7)) + ...
quad_form(p_9,psi_sck(:,:,7)) + t_sck(7) - 2*real(f_sck(:,7)'*p_sc) + u_sck(7) - v_sck(7);

X_sc8 = quad_form(p_sc,psi_sck(:,:,8)) + quad_form(p_c1,psi_sck(:,:,8)) + quad_form(p_c2,psi_sck(:,:,8)) + ...
quad_form(p_c3,psi_sck(:,:,8)) + quad_form(p_1,psi_sck(:,:,8)) + quad_form(p_2,psi_sck(:,:,8)) + ...
quad_form(p_3,psi_sck(:,:,8)) + quad_form(p_4,psi_sck(:,:,8)) + quad_form(p_5,psi_sck(:,:,8)) + ...
quad_form(p_6,psi_sck(:,:,8)) + quad_form(p_7,psi_sck(:,:,8)) + quad_form(p_8,psi_sck(:,:,8)) + ...
quad_form(p_9,psi_sck(:,:,8)) + t_sck(8) - 2*real(f_sck(:,8)'*p_sc) + u_sck(8) - v_sck(8);

X_sc9 = quad_form(p_sc,psi_sck(:,:,9)) + quad_form(p_c1,psi_sck(:,:,9)) + quad_form(p_c2,psi_sck(:,:,9)) + ...
quad_form(p_c3,psi_sck(:,:,9)) + quad_form(p_1,psi_sck(:,:,9)) + quad_form(p_2,psi_sck(:,:,9)) + ...
quad_form(p_3,psi_sck(:,:,9)) + quad_form(p_4,psi_sck(:,:,9)) + quad_form(p_5,psi_sck(:,:,9)) + ...
quad_form(p_6,psi_sck(:,:,9)) + quad_form(p_7,psi_sck(:,:,9)) + quad_form(p_8,psi_sck(:,:,9)) + ...
quad_form(p_9,psi_sck(:,:,9)) + t_sck(9) - 2*real(f_sck(:,9)'*p_sc) + u_sck(9) - v_sck(9);

X_sc10 = quad_form(p_sc,psi_sck(:,:,10)) + quad_form(p_c1,psi_sck(:,:,10)) + quad_form(p_c2,psi_sck(:,:,10)) + ...
quad_form(p_c3,psi_sck(:,:,10)) + quad_form(p_1,psi_sck(:,:,10)) + quad_form(p_2,psi_sck(:,:,10)) + ...
quad_form(p_3,psi_sck(:,:,10)) + quad_form(p_4,psi_sck(:,:,10)) + quad_form(p_5,psi_sck(:,:,10)) + ...
quad_form(p_6,psi_sck(:,:,10)) + quad_form(p_7,psi_sck(:,:,10)) + quad_form(p_8,psi_sck(:,:,10)) + ...
quad_form(p_9,psi_sck(:,:,10)) + t_sck(10) - 2*real(f_sck(:,10)'*p_sc) + u_sck(10) - v_sck(10);

X_sc11 = quad_form(p_sc,psi_sck(:,:,11)) + quad_form(p_c1,psi_sck(:,:,11)) + quad_form(p_c2,psi_sck(:,:,11)) + ...
quad_form(p_c3,psi_sck(:,:,11)) + quad_form(p_1,psi_sck(:,:,11)) + quad_form(p_2,psi_sck(:,:,11)) + ...
quad_form(p_3,psi_sck(:,:,11)) + quad_form(p_4,psi_sck(:,:,11)) + quad_form(p_5,psi_sck(:,:,11)) + ...
quad_form(p_6,psi_sck(:,:,11)) + quad_form(p_7,psi_sck(:,:,11)) + quad_form(p_8,psi_sck(:,:,11)) + ...
quad_form(p_9,psi_sck(:,:,11)) + t_sck(11) - 2*real(f_sck(:,11)'*p_sc) + u_sck(11) - v_sck(11);

X_sc12 = quad_form(p_sc,psi_sck(:,:,12)) + quad_form(p_c1,psi_sck(:,:,12)) + quad_form(p_c2,psi_sck(:,:,12)) + ...
quad_form(p_c3,psi_sck(:,:,12)) + quad_form(p_1,psi_sck(:,:,12)) + quad_form(p_2,psi_sck(:,:,12)) + ...
quad_form(p_3,psi_sck(:,:,12)) + quad_form(p_4,psi_sck(:,:,12)) + quad_form(p_5,psi_sck(:,:,12)) + ...
quad_form(p_6,psi_sck(:,:,12)) + quad_form(p_7,psi_sck(:,:,12)) + quad_form(p_8,psi_sck(:,:,12)) + ...
quad_form(p_9,psi_sck(:,:,12)) + t_sck(12) - 2*real(f_sck(:,12)'*p_sc) + u_sck(12) - v_sck(12);

X_sc13 = quad_form(p_sc,psi_sck(:,:,13)) + quad_form(p_c1,psi_sck(:,:,13)) + quad_form(p_c2,psi_sck(:,:,13)) + ...
quad_form(p_c3,psi_sck(:,:,13)) + quad_form(p_1,psi_sck(:,:,13)) + quad_form(p_2,psi_sck(:,:,13)) + ...
quad_form(p_3,psi_sck(:,:,13)) + quad_form(p_4,psi_sck(:,:,13)) + quad_form(p_5,psi_sck(:,:,13)) + ...
quad_form(p_6,psi_sck(:,:,13)) + quad_form(p_7,psi_sck(:,:,13)) + quad_form(p_8,psi_sck(:,:,13)) + ...
quad_form(p_9,psi_sck(:,:,13)) + t_sck(13) - 2*real(f_sck(:,13)'*p_sc) + u_sck(13) - v_sck(13);

X_sc14 = quad_form(p_sc,psi_sck(:,:,14)) + quad_form(p_c1,psi_sck(:,:,14)) + quad_form(p_c2,psi_sck(:,:,14)) + ...
quad_form(p_c3,psi_sck(:,:,14)) + quad_form(p_1,psi_sck(:,:,14)) + quad_form(p_2,psi_sck(:,:,14)) + ...
quad_form(p_3,psi_sck(:,:,14)) + quad_form(p_4,psi_sck(:,:,14)) + quad_form(p_5,psi_sck(:,:,14)) + ...
quad_form(p_6,psi_sck(:,:,14)) + quad_form(p_7,psi_sck(:,:,14)) + quad_form(p_8,psi_sck(:,:,14)) + ...
quad_form(p_9,psi_sck(:,:,14)) + t_sck(14) - 2*real(f_sck(:,14)'*p_sc) + u_sck(14) - v_sck(14);

X_sc15 = quad_form(p_sc,psi_sck(:,:,15)) + quad_form(p_c1,psi_sck(:,:,15)) + quad_form(p_c2,psi_sck(:,:,15)) + ...
quad_form(p_c3,psi_sck(:,:,15)) + quad_form(p_1,psi_sck(:,:,15)) + quad_form(p_2,psi_sck(:,:,15)) + ...
quad_form(p_3,psi_sck(:,:,15)) + quad_form(p_4,psi_sck(:,:,15)) + quad_form(p_5,psi_sck(:,:,15)) + ...
quad_form(p_6,psi_sck(:,:,15)) + quad_form(p_7,psi_sck(:,:,15)) + quad_form(p_8,psi_sck(:,:,15)) + ...
quad_form(p_9,psi_sck(:,:,15)) + t_sck(15) - 2*real(f_sck(:,15)'*p_sc) + u_sck(15) - v_sck(15);

X_sc16 = quad_form(p_sc,psi_sck(:,:,16)) + quad_form(p_c1,psi_sck(:,:,16)) + quad_form(p_c2,psi_sck(:,:,16)) + ...
quad_form(p_c3,psi_sck(:,:,16)) + quad_form(p_1,psi_sck(:,:,16)) + quad_form(p_2,psi_sck(:,:,16)) + ...
quad_form(p_3,psi_sck(:,:,16)) + quad_form(p_4,psi_sck(:,:,16)) + quad_form(p_5,psi_sck(:,:,16)) + ...
quad_form(p_6,psi_sck(:,:,16)) + quad_form(p_7,psi_sck(:,:,16)) + quad_form(p_8,psi_sck(:,:,16)) + ...
quad_form(p_9,psi_sck(:,:,16)) + t_sck(16) - 2*real(f_sck(:,16)'*p_sc) + u_sck(16) - v_sck(16);

X_sc17 = quad_form(p_sc,psi_sck(:,:,17)) + quad_form(p_c1,psi_sck(:,:,17)) + quad_form(p_c2,psi_sck(:,:,17)) + ...
quad_form(p_c3,psi_sck(:,:,17)) + quad_form(p_1,psi_sck(:,:,17)) + quad_form(p_2,psi_sck(:,:,17)) + ...
quad_form(p_3,psi_sck(:,:,17)) + quad_form(p_4,psi_sck(:,:,17)) + quad_form(p_5,psi_sck(:,:,17)) + ...
quad_form(p_6,psi_sck(:,:,17)) + quad_form(p_7,psi_sck(:,:,17)) + quad_form(p_8,psi_sck(:,:,17)) + ...
quad_form(p_9,psi_sck(:,:,17)) + t_sck(17) - 2*real(f_sck(:,17)'*p_sc) + u_sck(17) - v_sck(17);

X_sc18 = quad_form(p_sc,psi_sck(:,:,18)) + quad_form(p_c1,psi_sck(:,:,18)) + quad_form(p_c2,psi_sck(:,:,18)) + ...
quad_form(p_c3,psi_sck(:,:,18)) + quad_form(p_1,psi_sck(:,:,18)) + quad_form(p_2,psi_sck(:,:,18)) + ...
quad_form(p_3,psi_sck(:,:,18)) + quad_form(p_4,psi_sck(:,:,18)) + quad_form(p_5,psi_sck(:,:,18)) + ...
quad_form(p_6,psi_sck(:,:,18)) + quad_form(p_7,psi_sck(:,:,18)) + quad_form(p_8,psi_sck(:,:,18)) + ...
quad_form(p_9,psi_sck(:,:,18)) + t_sck(18) - 2*real(f_sck(:,18)'*p_sc) + u_sck(18) - v_sck(18);

%Objective Function
object_func = r_g;

%Optimisation
maximize(object_func)

%Constraints    
constraints(1) = SC_1 + C_1 + r_1 - r_g;
constraints(2) = SC_2 + C_2 + r_2 - r_g;
constraints(3) = SC_3 + C_3 + r_3 - r_g;
constraints(4) = SC_4 + C_4 + r_4 - r_g;
constraints(5) = SC_5 + C_5 + r_5 - r_g;
constraints(6) = SC_6 + C_6 + r_6 - r_g;
constraints(7) = SC_7 + C_7 + r_7 - r_g;
constraints(8) = SC_8 + C_8 + r_8 - r_g;
constraints(9) = SC_9 + C_9 + r_9 - r_g;
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
constraints(46) = 1 - X_sc1 - SC_1 - SC_2 - SC_3 - SC_4 - SC_5 - SC_6 - SC_7 - SC_8 - SC_9;
constraints(47) = 1 - X_sc2 - SC_1 - SC_2 - SC_3 - SC_4 - SC_5 - SC_6 - SC_7 - SC_8 - SC_9;
constraints(48) = 1 - X_sc3 - SC_1 - SC_2 - SC_3 - SC_4 - SC_5 - SC_6 - SC_7 - SC_8 - SC_9;
constraints(49) = 1 - X_sc4 - SC_1 - SC_2 - SC_3 - SC_4 - SC_5 - SC_6 - SC_7 - SC_8 - SC_9;
constraints(50) = 1 - X_sc5 - SC_1 - SC_2 - SC_3 - SC_4 - SC_5 - SC_6 - SC_7 - SC_8 - SC_9;
constraints(51) = 1 - X_sc6 - SC_1 - SC_2 - SC_3 - SC_4 - SC_5 - SC_6 - SC_7 - SC_8 - SC_9;
constraints(52) = 1 - X_sc7 - SC_1 - SC_2 - SC_3 - SC_4 - SC_5 - SC_6 - SC_7 - SC_8 - SC_9;
constraints(53) = 1 - X_sc8 - SC_1 - SC_2 - SC_3 - SC_4 - SC_5 - SC_6 - SC_7 - SC_8 - SC_9;
constraints(54) = 1 - X_sc9 - SC_1 - SC_2 - SC_3 - SC_4 - SC_5 - SC_6 - SC_7 - SC_8 - SC_9;
constraints(55) = 1 - X_sc10 - SC_1 - SC_2 - SC_3 - SC_4 - SC_5 - SC_6 - SC_7 - SC_8 - SC_9;
constraints(56) = 1 - X_sc11 - SC_1 - SC_2 - SC_3 - SC_4 - SC_5 - SC_6 - SC_7 - SC_8 - SC_9;
constraints(57) = 1 - X_sc12 - SC_1 - SC_2 - SC_3 - SC_4 - SC_5 - SC_6 - SC_7 - SC_8 - SC_9;
constraints(58) = 1 - X_sc13 - SC_1 - SC_2 - SC_3 - SC_4 - SC_5 - SC_6 - SC_7 - SC_8 - SC_9;
constraints(59) = 1 - X_sc14 - SC_1 - SC_2 - SC_3 - SC_4 - SC_5 - SC_6 - SC_7 - SC_8 - SC_9;
constraints(60) = 1 - X_sc15 - SC_1 - SC_2 - SC_3 - SC_4 - SC_5 - SC_6 - SC_7 - SC_8 - SC_9;
constraints(61) = 1 - X_sc16 - SC_1 - SC_2 - SC_3 - SC_4 - SC_5 - SC_6 - SC_7 - SC_8 - SC_9;
constraints(62) = 1 - X_sc17 - SC_1 - SC_2 - SC_3 - SC_4 - SC_5 - SC_6 - SC_7 - SC_8 - SC_9;
constraints(63) = 1 - X_sc18 - SC_1 - SC_2 - SC_3 - SC_4 - SC_5 - SC_6 - SC_7 - SC_8 - SC_9;
constraints(64) = C_1;
constraints(65) = C_2;
constraints(66) = C_3;
constraints(67) = C_4;
constraints(68) = C_5;
constraints(69) = C_6;
constraints(70) = C_7;
constraints(71) = C_8;
constraints(72) = C_9;
constraints(73) = SC_1;
constraints(74) = SC_2;
constraints(75) = SC_3;
constraints(76) = SC_4;
constraints(77) = SC_5;
constraints(78) = SC_6;
constraints(79) = SC_7;
constraints(80) = SC_8;
constraints(81) = SC_9;

constraints(82) = PAC - (p_sc'*D1*p_sc + p_c1'*D1*p_c1 + ...
p_c2'*D1*p_c2 + p_c3'*D1*p_c3 + p_1'*D1*p_1 + p_2'*D1*p_2 + p_3'*D1*p_3 + ...
p_4'*D1*p_4 + p_5'*D1*p_5 + p_6'*D1*p_6 + p_7'*D1*p_7 + p_8'*D1*p_8 + p_9'*D1*p_9);

constraints(83) = PAC - (p_sc'*D2*p_sc + p_c1'*D2*p_c1 + ...
p_c2'*D2*p_c2 + p_c3'*D2*p_c3 + p_1'*D2*p_1 + p_2'*D2*p_2 + p_3'*D2*p_3 + ...
p_4'*D2*p_4 + p_5'*D2*p_5 + p_6'*D2*p_6 + p_7'*D2*p_7 + p_8'*D2*p_8 + p_9'*D2*p_9);

constraints(84) = PAC - (p_sc'*D3*p_sc + p_c1'*D3*p_c1 + ...
p_c2'*D3*p_c2 + p_c3'*D3*p_c3 + p_1'*D3*p_1 + p_2'*D3*p_2 + p_3'*D3*p_3 +...
p_4'*D3*p_4 + p_5'*D3*p_5 + p_6'*D3*p_6 + p_7'*D3*p_7 + p_8'*D3*p_8 + p_9'*D3*p_9);

constraints(85) = PAC - (p_sc'*D4*p_sc + p_c1'*D4*p_c1 + ...
p_c2'*D4*p_c2 + p_c3'*D4*p_c3 + p_1'*D4*p_1 + p_2'*D4*p_2 + p_3'*D4*p_3 +...
p_4'*D4*p_4 + p_5'*D4*p_5 + p_6'*D4*p_6 + p_7'*D4*p_7 + p_8'*D4*p_8 + p_9'*D4*p_9);

constraints(86) = PAC - (p_sc'*D5*p_sc + p_c1'*D5*p_c1 + ...
p_c2'*D5*p_c2 + p_c3'*D5*p_c3 + p_1'*D5*p_1 + p_2'*D5*p_2 + p_3'*D5*p_3 +...
p_4'*D5*p_4 + p_5'*D5*p_5 + p_6'*D5*p_6 + p_7'*D5*p_7 + p_8'*D5*p_8 + p_9'*D5*p_9);

constraints(87) = PAC - (p_sc'*D6*p_sc + p_c1'*D6*p_c1 + ...
p_c2'*D6*p_c2 + p_c3'*D6*p_c3 + p_1'*D6*p_1 + p_2'*D6*p_2 + p_3'*D6*p_3 +...
p_4'*D6*p_4 + p_5'*D6*p_5 + p_6'*D6*p_6 + p_7'*D6*p_7 + p_8'*D6*p_8 + p_9'*D6*p_9);

constraints(88) = PAC - (p_sc'*D7*p_sc + p_c1'*D7*p_c1 + ...
p_c2'*D7*p_c2 + p_c3'*D7*p_c3 + p_1'*D7*p_1 + p_2'*D7*p_2 + p_3'*D7*p_3 +...
p_4'*D7*p_4 + p_5'*D7*p_5 + p_6'*D7*p_6 + p_7'*D7*p_7 + p_8'*D7*p_8 + p_9'*D7*p_9);

constraints(89) = PAC - (p_sc'*D8*p_sc + p_c1'*D8*p_c1 + ...
p_c2'*D8*p_c2 + p_c3'*D8*p_c3 + p_1'*D8*p_1 + p_2'*D8*p_2 + p_3'*D8*p_3 +...
p_4'*D8*p_4 + p_5'*D8*p_5 + p_6'*D8*p_6 + p_7'*D8*p_7 + p_8'*D8*p_8 + p_9'*D8*p_9);

constraints(90) = PAC - (p_sc'*D9*p_sc + p_c1'*D9*p_c1 + ...
p_c2'*D9*p_c2 + p_c3'*D9*p_c3 + p_1'*D9*p_1 + p_2'*D9*p_2 + p_3'*D9*p_3 +...
p_4'*D9*p_4 + p_5'*D9*p_5 + p_6'*D9*p_6 + p_7'*D9*p_7 + p_8'*D9*p_8 + p_9'*D9*p_9);

subject to
    constraints >= zeros(size(constraints))

cvx_end

%Result
rate = object_func;
