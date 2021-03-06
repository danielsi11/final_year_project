function MMF = NoRS_perfect_rate(Gm,H,PAC,tolerance)

[B,K,L] = size(H);
P = PAC*B*L; 

%Initialisation
power = P/(B*L);

%Precoder for Private Stream
offset = 0;
for i_cluster = 1:L
    for i_group = 1:B
        H_Gm = [];
        for i_user = 1:Gm(i_group,i_cluster)
            H_Gm = horzcat(H_Gm,H(:,i_user+offset,i_cluster));
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
    
    %Update equalisers and weights
    %Received Power
    %T_k Calculation
    T_k = zeros(1,K);
    for i_user = 1:K         
        T_k(i_user) = T_k(i_user) + 1;
        
        for i_cluster = 1:L
            for i_group = 1:B   
                T_k(i_user) = T_k(i_user) + abs(H(:,i_user,i_cluster)'*p_m(:,i_group,i_cluster))^2;
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
        
        g_k(i_user) = p_m(:,i_group,i_cluster)'*H(:,i_user,i_cluster)*inv(T_k(i_user));
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
        
        E_k(i_user) = inv(T_k(i_user))*(T_k(i_user)-abs(H(:,i_user,i_cluster)'*p_m(:,i_group,i_cluster))^2);
    end
    
    %MMSE Weights
    %W_k Calculation
    W_k = zeros(1,K);
    for i_user = 1:K
        W_k(i_user) = inv(E_k(i_user));
    end
    
    %Update P and Find MMF rate
    [rate, p_1, p_2, p_3, p_4, p_5, p_6, p_7, p_8, p_9] = NoRS_perfect_optimisation(H,PAC,g_k,W_k);
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
