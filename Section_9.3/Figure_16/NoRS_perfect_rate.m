function MMF = NoRS_perfect_rate(Gm,H,PAC,tolerance,G)

[B,K,L] = size(H);
P = PAC*B*L; 

%1st Stage Precoding
%Initialisation
rng('default');
rng(1);

for i_gateway = 1:L
    U(:,:,i_gateway) = 1/sqrt(2)*(randn(B,B)+1i*randn(B,B));
end

%Alternating Optimisation
loop=1;
sum_MSE_past=0;
count=0;
while (loop)
    
    %Update OBP matrix R
    for i_cluster = 1:L
        temp = eye(B);
        for i_gateway = 1:L
            temp = temp + G(:,:,i_cluster,i_gateway)*U(:,:,i_gateway)*U(:,:,i_gateway)'*G(:,:,i_cluster,i_gateway)';
        end

        R(:,:,i_cluster) = U(:,:,i_cluster)'*G(:,:,i_cluster,i_cluster)'*inv(temp);
    end

    %Update 1st stage precoder U
    for i_gateway = 1:L
        temp = zeros(B,B);
        for i_cluster = 1:L
            temp = temp + G(:,:,i_cluster,i_gateway)'*R(:,:,i_cluster)'*R(:,:,i_cluster)*G(:,:,i_cluster,i_gateway);
        end
        U(:,:,i_gateway) = inv(temp)*G(:,:,i_gateway,i_gateway)'*R(:,:,i_gateway)';
    end

    %Find MSE
    for i_cluster = 1:L
        temp = zeros(B,B);
        for i_gateway = 1:L
            temp = temp + G(:,:,i_cluster,i_gateway)*U(:,:,i_gateway)*U(:,:,i_gateway)'*G(:,:,i_cluster,i_gateway)';
        end

        MSE(i_cluster) = trace(R(:,:,i_cluster)*temp*R(:,:,i_cluster)' - R(:,:,i_cluster)*G(:,:,i_cluster,i_cluster)*U(:,:,i_cluster) - ...
        U(:,:,i_cluster)'*G(:,:,i_cluster,i_cluster)'*R(:,:,i_cluster)' + eye(B) + R(:,:,i_cluster)*R(:,:,i_cluster)');
    end

    sum_MSE = sum(MSE);
    
    if abs(sum_MSE - sum_MSE_past) <= tolerance
        loop = 0;
    else
        sum_MSE_past = sum_MSE;
        count = count+1;
    end
    
    if count >= 2000 %Exit loop if convergence has yet to be reached by this point
        break;
    end
    
end

%Effective Feeder Link Channel
for i_gateway = 1:L
    for i_cluster = 1:L
        G_eff(:,:,i_cluster,i_gateway) = R(:,:,i_cluster)*G(:,:,i_cluster,i_gateway)*U(:,:,i_gateway);
    end
end

%Overall Effective Channel
for i_user = 1:K
    for i_gateway = 1:L
        temp = zeros(1,B);
        for i_cluster = 1:L
            temp = temp + H(:,i_user,i_cluster)'*G_eff(:,:,i_cluster,i_gateway);
        end
        H_eff(:,i_user,i_gateway) = temp';
    end
end

%Overall Effective Noise Power
for i_user = 1:K    
    P_n_eff(i_user) = 1;
    for i_cluster = 1:L
        P_n_eff(i_user) = P_n_eff(i_user) + norm(H(:,i_user,i_cluster)'*R(:,:,i_cluster))^2;
    end    
end

%2nd Stage Precoding
%Initialisation
power = P/(B*L);

%Precoder for Private Stream
offset = 0;
for i_cluster = 1:L
    for i_group = 1:B
        H_Gm = [];
        for i_user = 1:Gm(i_group,i_cluster)
            H_Gm = horzcat(H_Gm,H_eff(:,i_user+offset,i_cluster));
        end
        [U_Gm,~,~] = svd(H_Gm);
        v_m(:,i_group,i_cluster) = U_Gm(:,1)*sqrt(power);
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
        T_k(i_user) = T_k(i_user) + P_n_eff(i_user);
        
        for i_cluster = 1:L
            for i_group = 1:B   
                T_k(i_user) = T_k(i_user) + abs(H_eff(:,i_user,i_cluster)'*v_m(:,i_group,i_cluster))^2;
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
        
        g_k(i_user) = v_m(:,i_group,i_cluster)'*H_eff(:,i_user,i_cluster)*inv(T_k(i_user));
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
        
        E_k(i_user) = inv(T_k(i_user))*(T_k(i_user)-abs(H_eff(:,i_user,i_cluster)'*v_m(:,i_group,i_cluster))^2);
    end

    %MMSE Weights
    %W_k Calculation
    W_k = zeros(1,K);
    for i_user = 1:K
        W_k(i_user) = inv(E_k(i_user));
    end
    
    %Update 2nd stage precoder V and Find MMF rate
    [rate, v_1, v_2, v_3, v_4, v_5, v_6, v_7, v_8, v_9] = NoRS_perfect_optimisation(H_eff,PAC,g_k,W_k,P_n_eff,G_eff,R);
    v_m(:,:,1) = horzcat(v_1, v_2, v_3);
    v_m(:,:,2) = horzcat(v_4, v_5, v_6);
    v_m(:,:,3) = horzcat(v_7, v_8, v_9);
    
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
