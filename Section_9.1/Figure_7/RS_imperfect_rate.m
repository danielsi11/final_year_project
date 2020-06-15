function MMF = RS_imperfect_rate(Gm,H_est,PAC,tolerance,S,P_e)

[N,K] = size(H_est); 
M = numel(Gm);
P = PAC*N;

%Channel Realisation Set
rng('default');
rng(1);

for i_sample = 1:S
    for i_user = 1:K
        H_error_sample(:,i_user) = sqrt(P_e)/sqrt(2)*(randn(N,1)+1i*randn(N,1));
    end    
    H(:,:,i_sample) = H_est + H_error_sample;
end

%Initialisation
power = P/(M + 1);

%Precoder for Global Common Stream
[U_sc,~,~] = svd(H_est);
p_sc = U_sc(:,1)*sqrt(power);

%Precoder for Private Stream
offset = 0;
for i_group = 1:M
    H_Gm = [];
    for i_user = 1:Gm(i_group)
        H_Gm = horzcat(H_Gm,H_est(:,i_user+offset));
    end
    [U_Gm,~,~] = svd(H_Gm);
    p_m(:,i_group) = U_Gm(:,1)*sqrt(power);
    offset = offset + Gm(i_group);
end

%Alternating Optimisation
loop=1;
rate_past=0;
count=0;
while (loop)
    
    t_k = zeros(1,K);
    t_sck = zeros(1,K);
    psi_k = zeros(N,N,K);
    psi_sck = zeros(N,N,K);
    f_k = zeros(N,K);
    f_sck = zeros(N,K);
    v_k = zeros(1,K);
    v_sck = zeros(1,K);
    u_k = zeros(1,K);
    u_sck = zeros(1,K);
    
    for i_sample = 1:S
        
        H_sample = H(:,:,i_sample);
    
        %Update equalisers and weights
        %Received Power
        %T_k Calculation
        T_k = zeros(1,K);
        for i_user = 1:K
            for i_group = 1:M
            T_k(i_user) = T_k(i_user) + abs(H_sample(:,i_user)'*p_m(:,i_group))^2;
            end
            T_k(i_user) = T_k(i_user) + 1;
        end

        %T_sck Calculation
        T_sck = zeros(1,K);
        for i_user = 1:K
            T_sck(i_user) = abs(H_sample(:,i_user)'*p_sc)^2 + T_k(i_user);
        end

        %MMSE Equalisers
        %g_k Calculation
        g_k = zeros(1,K);
        offset = 0;
        i_group = 1;
        for i_user = 1:K
            if i_user > Gm(i_group) + offset
                offset = offset + Gm(i_group);
                i_group = i_group + 1;
            end
            g_k(i_user) = p_m(:,i_group)'*H_sample(:,i_user)*inv(T_k(i_user));
        end

        %g_sck Calculation
        g_sck = zeros(1,K);
        for i_user = 1:K
            g_sck(i_user) = p_sc'*H_sample(:,i_user)*inv(T_sck(i_user));
        end

        %MMSE
        %E_k Calculation
        E_k = zeros(1,K);
        offset = 0;
        i_group = 1;
        for i_user = 1:K
            if i_user > Gm(i_group) + offset
                offset = offset + Gm(i_group);
                i_group = i_group + 1;
            end
            E_k(i_user) = inv(T_k(i_user))*(T_k(i_user)-abs(H_sample(:,i_user)'*p_m(:,i_group))^2);
        end

        %E_sck Calculation
        E_sck = zeros(1,K);
        for i_user = 1:K
            E_sck(i_user) = inv(T_sck(i_user))*(T_sck(i_user)-abs(H_sample(:,i_user)'*p_sc)^2);
        end

        %MMSE Weights
        %W_k Calculation
        W_k = zeros(1,K);
        for i_user = 1:K
            W_k(i_user) = inv(E_k(i_user));
        end
        
        %W_sck Calculation
        W_sck = zeros(1,K);
        for i_user = 1:K
            W_sck(i_user) = inv(E_sck(i_user));
        end
                
        %SAF Calculation
        for i_user = 1:K
            t_k(i_user) = t_k(i_user) + W_k(i_user)*abs(g_k(i_user))^2;
            t_sck(i_user) = t_sck(i_user) + W_sck(i_user)*abs(g_sck(i_user))^2;

            psi_k(:,:,i_user) = psi_k(:,:,i_user) + W_k(i_user)*abs(g_k(i_user))^2*H_sample(:,i_user)*H_sample(:,i_user)';
            psi_sck(:,:,i_user) = psi_sck(:,:,i_user) + W_sck(i_user)*abs(g_sck(i_user))^2*H_sample(:,i_user)*H_sample(:,i_user)';

            f_k(:,i_user) = f_k(:,i_user) + W_k(i_user)*H_sample(:,i_user)*g_k(i_user)';
            f_sck(:,i_user) = f_sck(:,i_user) + W_sck(i_user)*H_sample(:,i_user)*g_sck(i_user)';

            v_k(i_user) = v_k(i_user) + log2(W_k(i_user));
            v_sck(i_user) = v_sck(i_user) + log2(W_sck(i_user));

            u_k(i_user) = u_k(i_user) + W_k(i_user);
            u_sck(i_user) = u_sck(i_user) + W_sck(i_user);
        end    
        
    end
    
    t_k = t_k/S;
    t_sck = t_sck/S;
    psi_k = psi_k/S;
    psi_sck = psi_sck/S;
    f_k = f_k/S;
    f_sck = f_sck/S;
    v_k = v_k/S;
    v_sck = v_sck/S;
    u_k = u_k/S;
    u_sck = u_sck/S;

    %Update P and Find MMF rate
    [rate, p_sc, p_1, p_2, p_3, p_4, p_5, p_6, p_7] = RS_imperfect_optimisation(H_est,PAC,t_k,t_sck,psi_k,psi_sck,f_k,f_sck,v_k,v_sck,u_k,u_sck);
    p_m = horzcat(p_1, p_2, p_3, p_4, p_5, p_6, p_7);
    
    if abs(rate - rate_past) <= tolerance
        loop = 0;
    else
        rate_past = rate;
        count = count+1;
    end
    
    if count >= 2000 %Exit loop if convergence has yet to be reached by this point
        break;
    end
    
end

MMF = rate;
