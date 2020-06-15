function MMF = RS_perfect_rate(Gm,H,PAC,tolerance)

[N,K] = size(H);
M = numel(Gm);
P = PAC*N; 

%Initialisation
power = P/(M + 1);

%Precoder for Global Common Stream
[U_sc,~,~] = svd(H);
p_sc = U_sc(:,1)*sqrt(power);

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
    %Recieved Power
    %T_k Calculation
    T_k = zeros(1,K);
    for i_user = 1:K         
        for i_group = 1:M   
            T_k(i_user) = T_k(i_user) + abs(H(:,i_user)'*p_m(:,i_group))^2;
        end
        T_k(i_user) = T_k(i_user) + 1;
    end

    %T_sck Calculation
    T_sck = zeros(1,K);
    for i_user = 1:K
        T_sck(i_user) = abs(H(:,i_user)'*p_sc)^2 + T_k(i_user);
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
    
    %g_sck Calculation
    g_sck = zeros(1,K);
    for i_user = 1:K
        g_sck(i_user) = p_sc'*H(:,i_user)*inv(T_sck(i_user));
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
    
    %E_sck Calculation
    E_sck = zeros(1,K);
    for i_user = 1:K
        E_sck(i_user) = inv(T_sck(i_user))*(T_sck(i_user)-abs(H(:,i_user)'*p_sc)^2);
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

    %Update P and Find MMF rate
    [rate, p_sc, p_1, p_2, p_3, p_4, p_5, p_6, p_7, p_8, p_9] = RS_perfect_optimisation(H,PAC,g_k,g_sck,W_k,W_sck);
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
