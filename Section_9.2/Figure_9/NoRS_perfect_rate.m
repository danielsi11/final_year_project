function MMF = NoRS_perfect_rate(Gm,H,PAC,tolerance)

[N,K] = size(H);
M = numel(Gm);
P = PAC*N; 

%Initialisation
power = P/M;

%Precoder for Private Stream
offset = 0;
for i_group = 1:M
    H_Gm = [];
    for i_user = 1:Gm(i_group)
        H_Gm = horzcat(H_Gm,H(:,i_user+offset));
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
    
    %Update equalisers and weights
    %Received Power
    %T_k Calculation
    T_k = zeros(1,K);
    for i_user = 1:K         
        for i_group = 1:M   
            T_k(i_user) = T_k(i_user) + abs(H(:,i_user)'*p_m(:,i_group))^2;
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
        g_k(i_user) = p_m(:,i_group)'*H(:,i_user)*inv(T_k(i_user));
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
        E_k(i_user) = inv(T_k(i_user))*(T_k(i_user)-abs(H(:,i_user)'*p_m(:,i_group))^2);
    end

    %MMSE Weights
    %W_k Calculation
    W_k = zeros(1,K);
    for i_user = 1:K
        W_k(i_user) = inv(E_k(i_user));
    end
    
    %Update P and Find MMF rate
    [rate, p_1, p_2, p_3, p_4, p_5, p_6, p_7, p_8, p_9] = NoRS_perfect_optimisation(H,PAC,g_k,W_k);
    p_m = horzcat(p_1, p_2, p_3, p_4, p_5, p_6, p_7, p_8, p_9);
    
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
