function MMF = NoRS_imperfect_rate(Gm,H_est,PAC,tolerance,S,P_e)

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
power = P/M;

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
    psi_k = zeros(N,N,K);
    f_k = zeros(N,K);
    v_k = zeros(1,K);
    u_k = zeros(1,K);
    
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

        %MMSE Weights
        %W_k Calculation
        W_k = zeros(1,K);
        for i_user = 1:K
            W_k(i_user) = inv(E_k(i_user));
        end
                       
        %SAF Calculation
        for i_user = 1:K
            t_k(i_user) = t_k(i_user) + W_k(i_user)*abs(g_k(i_user))^2;
            psi_k(:,:,i_user) = psi_k(:,:,i_user) + W_k(i_user)*abs(g_k(i_user))^2*H_sample(:,i_user)*H_sample(:,i_user)';
            f_k(:,i_user) = f_k(:,i_user) + W_k(i_user)*H_sample(:,i_user)*g_k(i_user)';
            v_k(i_user) = v_k(i_user) + log2(W_k(i_user));
            u_k(i_user) = u_k(i_user) + W_k(i_user);
        end    
        
    end
    
    t_k = t_k/S;
    psi_k = psi_k/S;
    f_k = f_k/S;
    v_k = v_k/S;
    u_k = u_k/S;

    %Update P and Find MMF rate
    [rate, p_1, p_2, p_3, p_4, p_5, p_6, p_7] = NoRS_imperfect_optimisation(H_est,PAC,t_k,psi_k,f_k,v_k,u_k);
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
