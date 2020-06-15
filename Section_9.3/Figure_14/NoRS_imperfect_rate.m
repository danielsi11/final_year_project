function MMF = NoRS_imperfect_rate(Gm,H_est,PAC,tolerance,S,P_e)

[B,K,L] = size(H_est);
P = PAC*B*L;

%Channel Realisation Set
rng('default');
rng(1);

for i_sample = 1:S
    for i_user = 1:K
        for i_cluster = 1:L
            H_error_sample(:,i_user,i_cluster) = sqrt(P_e)/sqrt(2)*(randn(B,1)+1i*randn(B,1));
        end    
    end    
    H(:,:,:,i_sample) = H_est + H_error_sample;
end

%Initialisation
power = P/(B*L);

%Precoder for Private Stream
offset = 0;
for i_cluster = 1:L
    for i_group = 1:B
        H_Gm = [];
        for i_user = 1:Gm(i_group,i_cluster)
            H_Gm = horzcat(H_Gm,H_est(:,i_user+offset,i_cluster));
        end
        [U_Gm,~,~] = svd(H_Gm);
        p_m(:,i_group,i_cluster) = U_Gm(:,1)*sqrt(power);
        offset = offset + Gm(i_group,i_cluster);
    end    
end

%Alternating Optimisation
loop=1;
rate_past=0;
count=0;
while (loop)
    
    t_k = zeros(1,K);
    psi_k = zeros(B,B,K,L);
    f_k = zeros(B,K);
    v_k = zeros(1,K);
    u_k = zeros(1,K);
    
    for i_sample = 1:S
        
        H_sample = H(:,:,:,i_sample);
    
        %Update equalisers and weights
        %Received Power
        %T_k Calculation
        T_k = zeros(1,K);
        for i_user = 1:K         
            T_k(i_user) = T_k(i_user) + 1;

            for i_cluster = 1:L
                for i_group = 1:B   
                    T_k(i_user) = T_k(i_user) + abs(H_sample(:,i_user,i_cluster)'*p_m(:,i_group,i_cluster))^2;
                end
            end
        end    

        %MMSE Equalisers
        %g_k Calculation
        g_k = zeros(1,K);
        offset_group = 0;
        offset_cluster = 0;
        i_group = 1;
        i_cluster = 1;
        for i_user = 1:K

            if(i_user > sum(Gm(:,i_cluster)) + offset_cluster)
                offset_cluster = offset_cluster + sum(Gm(:,i_cluster));
                i_cluster = i_cluster + 1;
                offset_group = 0;
                i_group = 1;
            end    

            if(i_user > Gm(i_group,i_cluster) + offset_group + offset_cluster)
                offset_group = offset_group + Gm(i_group);
                i_group = i_group + 1;
            end

            g_k(i_user) = p_m(:,i_group,i_cluster)'*H_sample(:,i_user,i_cluster)*inv(T_k(i_user));
        end

        %MMSE
        %E_k Calculation
        E_k = zeros(1,K);
        offset_group = 0;
        offset_cluster = 0;
        i_group = 1;
        i_cluster = 1;
        for i_user = 1:K

            if(i_user > sum(Gm(:,i_cluster)) + offset_cluster)
                offset_cluster = offset_cluster + sum(Gm(:,i_cluster));
                i_cluster = i_cluster + 1;
                offset_group = 0;
                i_group = 1;
            end    

            if(i_user > Gm(i_group,i_cluster) + offset_group + offset_cluster)
                offset_group = offset_group + Gm(i_group);
                i_group = i_group + 1;
            end

            E_k(i_user) = inv(T_k(i_user))*(T_k(i_user)-abs(H_sample(:,i_user,i_cluster)'*p_m(:,i_group,i_cluster))^2);
        end

        %MMSE Weights
        %W_k Calculation
        W_k = zeros(1,K);
        for i_user = 1:K
            W_k(i_user) = inv(E_k(i_user));
        end
                       
        %SAF Calculation
        offset = 0;
        i_cluster = 1;
        for i_user = 1:K
            
            if i_user > sum(Gm(:,i_cluster)) + offset
                offset = offset + sum(Gm(:,i_cluster));
                i_cluster = i_cluster + 1;
            end
            
            t_k(i_user) = t_k(i_user) + W_k(i_user)*abs(g_k(i_user))^2;
            f_k(:,i_user) = f_k(:,i_user) + W_k(i_user)*H_sample(:,i_user,i_cluster)*g_k(i_user)';
            v_k(i_user) = v_k(i_user) + log2(W_k(i_user));
            u_k(i_user) = u_k(i_user) + W_k(i_user);
        end
        
        for i_cluster = 1:L
            for i_user = 1:K
                psi_k(:,:,i_user,i_cluster) = psi_k(:,:,i_user,i_cluster) + W_k(i_user)*abs(g_k(i_user))^2*H_sample(:,i_user,i_cluster)*H_sample(:,i_user,i_cluster)';
            end
        end    
        
    end
    
    t_k = t_k/S;
    psi_k = psi_k/S;
    f_k = f_k/S;
    v_k = v_k/S;
    u_k = u_k/S;

    %Update P and Find MMF rate
    [rate, p_1, p_2, p_3, p_4, p_5, p_6, p_7, p_8, p_9] = NoRS_imperfect_optimisation(H_est,PAC,t_k,psi_k,f_k,v_k,u_k);
    p_m(:,:,1) = horzcat(p_1, p_2, p_3);
    p_m(:,:,2) = horzcat(p_4, p_5, p_6);
    p_m(:,:,3) = horzcat(p_7, p_8, p_9);
    
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
