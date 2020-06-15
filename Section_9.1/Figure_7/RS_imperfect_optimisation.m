function [rate,p_sc,p_1,p_2,p_3,p_4,p_5,p_6,p_7] = RS_imperfect_optimisation(H_est,PAC,t_k,t_sck,psi_k,psi_sck,f_k,f_sck,v_k,v_sck,u_k,u_sck)
    
[N,K] = size(H_est); 

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
variable p_sc(N,1) complex
variable SC_1 
variable SC_2 
variable SC_3
variable SC_4 
variable SC_5 
variable SC_6
variable SC_7
variable r_1 
variable r_2 
variable r_3
variable r_4 
variable r_5 
variable r_6
variable r_7
variable r_g

expression constraints(1,49);

%Rate-WMMSE Relationship
%Private Stream
X_1 = quad_form(p_1,psi_k(:,:,1)) + quad_form(p_2,psi_k(:,:,1)) + quad_form(p_3,psi_k(:,:,1)) + ...
quad_form(p_4,psi_k(:,:,1)) + quad_form(p_5,psi_k(:,:,1)) + quad_form(p_6,psi_k(:,:,1)) + ...
quad_form(p_7,psi_k(:,:,1)) + t_k(1) - 2*real(f_k(:,1)'*p_1) + u_k(1) - v_k(1);

X_2 = quad_form(p_1,psi_k(:,:,2)) + quad_form(p_2,psi_k(:,:,2)) + quad_form(p_3,psi_k(:,:,2)) + ...
quad_form(p_4,psi_k(:,:,2)) + quad_form(p_5,psi_k(:,:,2)) + quad_form(p_6,psi_k(:,:,2)) + ...
quad_form(p_7,psi_k(:,:,2)) + t_k(2) - 2*real(f_k(:,2)'*p_1) + u_k(2) - v_k(2);

X_3 = quad_form(p_1,psi_k(:,:,3)) + quad_form(p_2,psi_k(:,:,3)) + quad_form(p_3,psi_k(:,:,3)) + ...
quad_form(p_4,psi_k(:,:,3)) + quad_form(p_5,psi_k(:,:,3)) + quad_form(p_6,psi_k(:,:,3)) + ...
quad_form(p_7,psi_k(:,:,3)) + t_k(3) - 2*real(f_k(:,3)'*p_2) + u_k(3) - v_k(3);

X_4 = quad_form(p_1,psi_k(:,:,4)) + quad_form(p_2,psi_k(:,:,4)) + quad_form(p_3,psi_k(:,:,4)) + ...
quad_form(p_4,psi_k(:,:,4)) + quad_form(p_5,psi_k(:,:,4)) + quad_form(p_6,psi_k(:,:,4)) + ...
quad_form(p_7,psi_k(:,:,4)) + t_k(4) - 2*real(f_k(:,4)'*p_2) + u_k(4) - v_k(4);

X_5 = quad_form(p_1,psi_k(:,:,5)) + quad_form(p_2,psi_k(:,:,5)) + quad_form(p_3,psi_k(:,:,5)) + ...
quad_form(p_4,psi_k(:,:,5)) + quad_form(p_5,psi_k(:,:,5)) + quad_form(p_6,psi_k(:,:,5)) + ...
quad_form(p_7,psi_k(:,:,5)) + t_k(5) - 2*real(f_k(:,5)'*p_3) + u_k(5) - v_k(5);

X_6 = quad_form(p_1,psi_k(:,:,6)) + quad_form(p_2,psi_k(:,:,6)) + quad_form(p_3,psi_k(:,:,6)) + ...
quad_form(p_4,psi_k(:,:,6)) + quad_form(p_5,psi_k(:,:,6)) + quad_form(p_6,psi_k(:,:,6)) + ...
quad_form(p_7,psi_k(:,:,6)) + t_k(6) - 2*real(f_k(:,6)'*p_3) + u_k(6) - v_k(6);

X_7 = quad_form(p_1,psi_k(:,:,7)) + quad_form(p_2,psi_k(:,:,7)) + quad_form(p_3,psi_k(:,:,7)) + ...
quad_form(p_4,psi_k(:,:,7)) + quad_form(p_5,psi_k(:,:,7)) + quad_form(p_6,psi_k(:,:,7)) + ...
quad_form(p_7,psi_k(:,:,7)) + t_k(7) - 2*real(f_k(:,7)'*p_4) + u_k(7) - v_k(7);

X_8 = quad_form(p_1,psi_k(:,:,8)) + quad_form(p_2,psi_k(:,:,8)) + quad_form(p_3,psi_k(:,:,8)) + ...
quad_form(p_4,psi_k(:,:,8)) + quad_form(p_5,psi_k(:,:,8)) + quad_form(p_6,psi_k(:,:,8)) + ...
quad_form(p_7,psi_k(:,:,8)) + t_k(8) - 2*real(f_k(:,8)'*p_4) + u_k(8) - v_k(8);

X_9 = quad_form(p_1,psi_k(:,:,9)) + quad_form(p_2,psi_k(:,:,9)) + quad_form(p_3,psi_k(:,:,9)) + ...
quad_form(p_4,psi_k(:,:,9)) + quad_form(p_5,psi_k(:,:,9)) + quad_form(p_6,psi_k(:,:,9)) + ...
quad_form(p_7,psi_k(:,:,9)) + t_k(9) - 2*real(f_k(:,9)'*p_5) + u_k(9) - v_k(9);

X_10 = quad_form(p_1,psi_k(:,:,10)) + quad_form(p_2,psi_k(:,:,10)) + quad_form(p_3,psi_k(:,:,10)) + ...
quad_form(p_4,psi_k(:,:,10)) + quad_form(p_5,psi_k(:,:,10)) + quad_form(p_6,psi_k(:,:,10)) + ...
quad_form(p_7,psi_k(:,:,10)) + t_k(10) - 2*real(f_k(:,10)'*p_5) + u_k(10) - v_k(10);

X_11 = quad_form(p_1,psi_k(:,:,11)) + quad_form(p_2,psi_k(:,:,11)) + quad_form(p_3,psi_k(:,:,11)) + ...
quad_form(p_4,psi_k(:,:,11)) + quad_form(p_5,psi_k(:,:,11)) + quad_form(p_6,psi_k(:,:,11)) + ...
quad_form(p_7,psi_k(:,:,11)) + t_k(11) - 2*real(f_k(:,11)'*p_6) + u_k(11) - v_k(11);

X_12 = quad_form(p_1,psi_k(:,:,12)) + quad_form(p_2,psi_k(:,:,12)) + quad_form(p_3,psi_k(:,:,12)) + ...
quad_form(p_4,psi_k(:,:,12)) + quad_form(p_5,psi_k(:,:,12)) + quad_form(p_6,psi_k(:,:,12)) + ...
quad_form(p_7,psi_k(:,:,12)) + t_k(12) - 2*real(f_k(:,12)'*p_6) + u_k(12) - v_k(12);

X_13 = quad_form(p_1,psi_k(:,:,13)) + quad_form(p_2,psi_k(:,:,13)) + quad_form(p_3,psi_k(:,:,13)) + ...
quad_form(p_4,psi_k(:,:,13)) + quad_form(p_5,psi_k(:,:,13)) + quad_form(p_6,psi_k(:,:,13)) + ...
quad_form(p_7,psi_k(:,:,13)) + t_k(13) - 2*real(f_k(:,13)'*p_7) + u_k(13) - v_k(13);

X_14 = quad_form(p_1,psi_k(:,:,14)) + quad_form(p_2,psi_k(:,:,14)) + quad_form(p_3,psi_k(:,:,14)) + ...
quad_form(p_4,psi_k(:,:,14)) + quad_form(p_5,psi_k(:,:,14)) + quad_form(p_6,psi_k(:,:,14)) + ...
quad_form(p_7,psi_k(:,:,14)) + t_k(14) - 2*real(f_k(:,14)'*p_7) + u_k(14) - v_k(14);

%Global Common Stream
X_sc1 = quad_form(p_sc,psi_sck(:,:,1)) + quad_form(p_1,psi_sck(:,:,1)) + quad_form(p_2,psi_sck(:,:,1)) + quad_form(p_3,psi_sck(:,:,1)) + ...
quad_form(p_4,psi_sck(:,:,1)) + quad_form(p_5,psi_sck(:,:,1)) + quad_form(p_6,psi_sck(:,:,1)) + ...
quad_form(p_7,psi_sck(:,:,1)) + t_sck(1) - 2*real(f_sck(:,1)'*p_sc) + u_sck(1) - v_sck(1);

X_sc2 = quad_form(p_sc,psi_sck(:,:,2)) + quad_form(p_1,psi_sck(:,:,2)) + quad_form(p_2,psi_sck(:,:,2)) + quad_form(p_3,psi_sck(:,:,2)) + ...
quad_form(p_4,psi_sck(:,:,2)) + quad_form(p_5,psi_sck(:,:,2)) + quad_form(p_6,psi_sck(:,:,2)) + ...
quad_form(p_7,psi_sck(:,:,2)) + t_sck(2) - 2*real(f_sck(:,2)'*p_sc) + u_sck(2) - v_sck(2);

X_sc3 = quad_form(p_sc,psi_sck(:,:,3)) + quad_form(p_1,psi_sck(:,:,3)) + quad_form(p_2,psi_sck(:,:,3)) + quad_form(p_3,psi_sck(:,:,3)) + ...
quad_form(p_4,psi_sck(:,:,3)) + quad_form(p_5,psi_sck(:,:,3)) + quad_form(p_6,psi_sck(:,:,3)) + ...
quad_form(p_7,psi_sck(:,:,3)) + t_sck(3) - 2*real(f_sck(:,3)'*p_sc) + u_sck(3) - v_sck(3);

X_sc4 = quad_form(p_sc,psi_sck(:,:,4)) + quad_form(p_1,psi_sck(:,:,4)) + quad_form(p_2,psi_sck(:,:,4)) + quad_form(p_3,psi_sck(:,:,4)) + ...
quad_form(p_4,psi_sck(:,:,4)) + quad_form(p_5,psi_sck(:,:,4)) + quad_form(p_6,psi_sck(:,:,4)) + ...
quad_form(p_7,psi_sck(:,:,4)) + t_sck(4) - 2*real(f_sck(:,4)'*p_sc) + u_sck(4) - v_sck(4);

X_sc5 = quad_form(p_sc,psi_sck(:,:,5)) + quad_form(p_1,psi_sck(:,:,5)) + quad_form(p_2,psi_sck(:,:,5)) + quad_form(p_3,psi_sck(:,:,5)) + ...
quad_form(p_4,psi_sck(:,:,5)) + quad_form(p_5,psi_sck(:,:,5)) + quad_form(p_6,psi_sck(:,:,5)) + ...
quad_form(p_7,psi_sck(:,:,5)) + t_sck(5) - 2*real(f_sck(:,5)'*p_sc) + u_sck(5) - v_sck(5);

X_sc6 = quad_form(p_sc,psi_sck(:,:,6)) + quad_form(p_1,psi_sck(:,:,6)) + quad_form(p_2,psi_sck(:,:,6)) + quad_form(p_3,psi_sck(:,:,6)) + ...
quad_form(p_4,psi_sck(:,:,6)) + quad_form(p_5,psi_sck(:,:,6)) + quad_form(p_6,psi_sck(:,:,6)) + ...
quad_form(p_7,psi_sck(:,:,6)) + t_sck(6) - 2*real(f_sck(:,6)'*p_sc) + u_sck(6) - v_sck(6);

X_sc7 = quad_form(p_sc,psi_sck(:,:,7)) + quad_form(p_1,psi_sck(:,:,7)) + quad_form(p_2,psi_sck(:,:,7)) + quad_form(p_3,psi_sck(:,:,7)) + ...
quad_form(p_4,psi_sck(:,:,7)) + quad_form(p_5,psi_sck(:,:,7)) + quad_form(p_6,psi_sck(:,:,7)) + ...
quad_form(p_7,psi_sck(:,:,7)) + t_sck(7) - 2*real(f_sck(:,7)'*p_sc) + u_sck(7) - v_sck(7);

X_sc8 = quad_form(p_sc,psi_sck(:,:,8)) + quad_form(p_1,psi_sck(:,:,8)) + quad_form(p_2,psi_sck(:,:,8)) + quad_form(p_3,psi_sck(:,:,8)) + ...
quad_form(p_4,psi_sck(:,:,8)) + quad_form(p_5,psi_sck(:,:,8)) + quad_form(p_6,psi_sck(:,:,8)) + ...
quad_form(p_7,psi_sck(:,:,8)) + t_sck(8) - 2*real(f_sck(:,8)'*p_sc) + u_sck(8) - v_sck(8);

X_sc9 = quad_form(p_sc,psi_sck(:,:,9)) + quad_form(p_1,psi_sck(:,:,9)) + quad_form(p_2,psi_sck(:,:,9)) + quad_form(p_3,psi_sck(:,:,9)) + ...
quad_form(p_4,psi_sck(:,:,9)) + quad_form(p_5,psi_sck(:,:,9)) + quad_form(p_6,psi_sck(:,:,9)) + ...
quad_form(p_7,psi_sck(:,:,9)) + t_sck(9) - 2*real(f_sck(:,9)'*p_sc) + u_sck(9) - v_sck(9);

X_sc10 = quad_form(p_sc,psi_sck(:,:,10)) + quad_form(p_1,psi_sck(:,:,10)) + quad_form(p_2,psi_sck(:,:,10)) + quad_form(p_3,psi_sck(:,:,10)) + ...
quad_form(p_4,psi_sck(:,:,10)) + quad_form(p_5,psi_sck(:,:,10)) + quad_form(p_6,psi_sck(:,:,10)) + ...
quad_form(p_7,psi_sck(:,:,10)) + t_sck(10) - 2*real(f_sck(:,10)'*p_sc) + u_sck(10) - v_sck(10);

X_sc11 = quad_form(p_sc,psi_sck(:,:,11)) + quad_form(p_1,psi_sck(:,:,11)) + quad_form(p_2,psi_sck(:,:,11)) + quad_form(p_3,psi_sck(:,:,11)) + ...
quad_form(p_4,psi_sck(:,:,11)) + quad_form(p_5,psi_sck(:,:,11)) + quad_form(p_6,psi_sck(:,:,11)) + ...
quad_form(p_7,psi_sck(:,:,11)) + t_sck(11) - 2*real(f_sck(:,11)'*p_sc) + u_sck(11) - v_sck(11);

X_sc12 = quad_form(p_sc,psi_sck(:,:,12)) + quad_form(p_1,psi_sck(:,:,12)) + quad_form(p_2,psi_sck(:,:,12)) + quad_form(p_3,psi_sck(:,:,12)) + ...
quad_form(p_4,psi_sck(:,:,12)) + quad_form(p_5,psi_sck(:,:,12)) + quad_form(p_6,psi_sck(:,:,12)) + ...
quad_form(p_7,psi_sck(:,:,12)) + t_sck(12) - 2*real(f_sck(:,12)'*p_sc) + u_sck(12) - v_sck(12);

X_sc13 = quad_form(p_sc,psi_sck(:,:,13)) + quad_form(p_1,psi_sck(:,:,13)) + quad_form(p_2,psi_sck(:,:,13)) + quad_form(p_3,psi_sck(:,:,13)) + ...
quad_form(p_4,psi_sck(:,:,13)) + quad_form(p_5,psi_sck(:,:,13)) + quad_form(p_6,psi_sck(:,:,13)) + ...
quad_form(p_7,psi_sck(:,:,13)) + t_sck(13) - 2*real(f_sck(:,13)'*p_sc) + u_sck(13) - v_sck(13);

X_sc14 = quad_form(p_sc,psi_sck(:,:,14)) + quad_form(p_1,psi_sck(:,:,14)) + quad_form(p_2,psi_sck(:,:,14)) + quad_form(p_3,psi_sck(:,:,14)) + ...
quad_form(p_4,psi_sck(:,:,14)) + quad_form(p_5,psi_sck(:,:,14)) + quad_form(p_6,psi_sck(:,:,14)) + ...
quad_form(p_7,psi_sck(:,:,14)) + t_sck(14) - 2*real(f_sck(:,14)'*p_sc) + u_sck(14) - v_sck(14);

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
constraints(22) = 1 - X_sc1 - SC_1 - SC_2 - SC_3 - SC_4 - SC_5 - SC_6 - SC_7;
constraints(23) = 1 - X_sc2 - SC_1 - SC_2 - SC_3 - SC_4 - SC_5 - SC_6 - SC_7;
constraints(24) = 1 - X_sc3 - SC_1 - SC_2 - SC_3 - SC_4 - SC_5 - SC_6 - SC_7;
constraints(25) = 1 - X_sc4 - SC_1 - SC_2 - SC_3 - SC_4 - SC_5 - SC_6 - SC_7;
constraints(26) = 1 - X_sc5 - SC_1 - SC_2 - SC_3 - SC_4 - SC_5 - SC_6 - SC_7;
constraints(27) = 1 - X_sc6 - SC_1 - SC_2 - SC_3 - SC_4 - SC_5 - SC_6 - SC_7;
constraints(28) = 1 - X_sc7 - SC_1 - SC_2 - SC_3 - SC_4 - SC_5 - SC_6 - SC_7;
constraints(29) = 1 - X_sc8 - SC_1 - SC_2 - SC_3 - SC_4 - SC_5 - SC_6 - SC_7;
constraints(30) = 1 - X_sc9 - SC_1 - SC_2 - SC_3 - SC_4 - SC_5 - SC_6 - SC_7;
constraints(31) = 1 - X_sc10 - SC_1 - SC_2 - SC_3 - SC_4 - SC_5 - SC_6 - SC_7;
constraints(32) = 1 - X_sc11 - SC_1 - SC_2 - SC_3 - SC_4 - SC_5 - SC_6 - SC_7;
constraints(33) = 1 - X_sc12 - SC_1 - SC_2 - SC_3 - SC_4 - SC_5 - SC_6 - SC_7;
constraints(34) = 1 - X_sc13 - SC_1 - SC_2 - SC_3 - SC_4 - SC_5 - SC_6 - SC_7;
constraints(35) = 1 - X_sc14 - SC_1 - SC_2 - SC_3 - SC_4 - SC_5 - SC_6 - SC_7;
constraints(36) = SC_1;
constraints(37) = SC_2;
constraints(38) = SC_3;
constraints(39) = SC_4;
constraints(40) = SC_5;
constraints(41) = SC_6;
constraints(42) = SC_7;
constraints(43) = PAC - (p_sc'*D1*p_sc + p_1'*D1*p_1 + p_2'*D1*p_2 + p_3'*D1*p_3 + ...
p_4'*D1*p_4 + p_5'*D1*p_5 + p_6'*D1*p_6 + p_7'*D1*p_7);
constraints(44) = PAC - (p_sc'*D2*p_sc + p_1'*D2*p_1 + p_2'*D2*p_2 + p_3'*D2*p_3 + ...
p_4'*D2*p_4 + p_5'*D2*p_5 + p_6'*D2*p_6 + p_7'*D2*p_7);
constraints(45) = PAC - (p_sc'*D3*p_sc + p_1'*D3*p_1 + p_2'*D3*p_2 + p_3'*D3*p_3 + ...
p_4'*D3*p_4 + p_5'*D3*p_5 + p_6'*D3*p_6 + p_7'*D3*p_7);
constraints(46) = PAC - (p_sc'*D4*p_sc + p_1'*D4*p_1 + p_2'*D4*p_2 + p_3'*D4*p_3 + ...
p_4'*D4*p_4 + p_5'*D4*p_5 + p_6'*D4*p_6 + p_7'*D4*p_7);
constraints(47) = PAC - (p_sc'*D5*p_sc + p_1'*D5*p_1 + p_2'*D5*p_2 + p_3'*D5*p_3 + ...
p_4'*D5*p_4 + p_5'*D5*p_5 + p_6'*D5*p_6 + p_7'*D5*p_7);
constraints(48) = PAC - (p_sc'*D6*p_sc + p_1'*D6*p_1 + p_2'*D6*p_2 + p_3'*D6*p_3 + ...
p_4'*D6*p_4 + p_5'*D6*p_5 + p_6'*D6*p_6 + p_7'*D6*p_7);
constraints(49) = PAC - (p_sc'*D7*p_sc + p_1'*D7*p_1 + p_2'*D7*p_2 + p_3'*D7*p_3 + ...
p_4'*D7*p_4 + p_5'*D7*p_5 + p_6'*D7*p_6 + p_7'*D7*p_7);

subject to
    constraints >= zeros(size(constraints))

cvx_end

%Result
rate = object_func;
