clc; clear all; clf

%Number of Transmitter Antennas
N = 9; 

%Number of Groups
M = 9;

%Number of Users per Group
rho = 2;

%Total Number of Users
K = rho*M;

%Number of Group Clusters
L = 3;

%Number of Users in Each Group
Gm = repmat(rho,M/L,L);

%Range of Per Antenna Power
PACs = 20:20:120; 

%Number of Random Channel
N_channel = 100;

%Accuracy of Convergence 
tolerance = 10^-4;

%Satellite System Parameter
Fc = 20*10^9;
height = 35786*10^3;
B_w = 500*10^6;
theta_3dB = 0.4*pi/180;
G_max = 10^(52/10);
G_R = 10^(41.7/10);
T_sys = 517;
rain_mu = -3.125;
rain_sigma = 1.591;
c = 3*10^8;
wavelength = c/Fc;
Boltz = 1.380649*10^-23;

%Satellite Location
satellite = [0;0;height];

%Beam Location
beam_radius = height*tan(theta_3dB);
theta = (0:1)*2*pi/6;
beam_spacing = ((beam_radius*3)/2)/cos(pi/6);

beam_temp = [beam_spacing*cos(theta); beam_spacing*sin(theta); zeros(1,2)];
beam_temp = horzcat(zeros(3,1),beam_temp);

shift0 = 1/2*(beam_temp(:,1)+beam_temp(:,3));
beam_temp = beam_temp - shift0;

shift1 = beam_temp(:,2) - beam_temp(:,3);
shift2 = beam_temp(:,1) - beam_temp(:,2);
shift3 = beam_temp(:,3) - beam_temp(:,1);

beam = horzcat(beam_temp + shift1, beam_temp + shift2, beam_temp + shift3);

%Simulate MMF rate
for i_channel = 1:N_channel
    
    for i_PAC = 1:length(PACs)

        rng('default');
        rng(i_channel+1);

        PAC = PACs(i_PAC);

        P = PAC*N; %Total Transmit Power
        
        offset = 0;
        i_group = 1;
        for i_user = 1:K

            %User Distribution
            if(i_user > Gm(i_group) + offset)
                offset = offset + Gm(i_group);
                i_group = i_group + 1;
            end    
            theta_rng = 2*pi*rand;
            user(:,i_user) = beam(:,i_group) + beam_radius*rand*[cos(theta_rng); sin(theta_rng); 0];
                       
            %User Link Channel Model
            X_k_dB(:,i_user) = exp(rain_mu + rain_sigma*randn*ones(N,1));
            X_k(:,i_user) = 10.^(X_k_dB(:,i_user)/20);

            phi_k(:,i_user) = 2*pi*rand*ones(N,1);

            q(:,i_user) = X_k(:,i_user).^(-0.5).*exp(-1i*phi_k(:,i_user));            

            d_k(i_user) = norm(satellite - user(:,i_user));

            user_sat = user(:,i_user) - satellite;

            for i_beam = 1:N
                beam_sat = beam(:,i_beam) - satellite;
                theta_nk(i_beam,i_user) = acos((beam_sat'*user_sat)/(norm(beam_sat)*norm(user_sat)));
            end 

            u_nk(:,i_user) = 2.07123*sin(theta_nk(:,i_user))/sin(theta_3dB);

            G_nk(:,i_user) = G_max*(besselj(1,u_nk(:,i_user))./(2*u_nk(:,i_user)) + ...
            36*besselj(3,u_nk(:,i_user))./(u_nk(:,i_user).^3)).^2;

            b(:,i_user) = sqrt(G_R*G_nk(:,i_user))./(4*pi*(d_k(i_user)/wavelength)*sqrt(Boltz*T_sys*B_w));

            H(:,i_user) = b(:,i_user).*q(:,i_user);

        end 

        %Perfect CSIT
        MMF_RS(i_PAC) = RS_perfect_rate(Gm,H,PAC,tolerance);
        MMF_NoRS(i_PAC) = NoRS_perfect_rate(Gm,H,PAC,tolerance);
        MMF_RS_cluster(i_PAC) = RS_cluster_perfect_rate(Gm,H,PAC,tolerance);
        MMF_HRS(i_PAC) = HRS_perfect_rate(Gm,H,PAC,tolerance);
        MMF_GRS(i_PAC) = GRS_perfect_rate(Gm,H,PAC,tolerance);

    end

    MMF_RS_set(i_channel,:) = MMF_RS;
    MMF_NoRS_set(i_channel,:) = MMF_NoRS;
    MMF_RS_cluster_set(i_channel,:) = MMF_RS_cluster;
    MMF_HRS_set(i_channel,:) = MMF_HRS;
    MMF_GRS_set(i_channel,:) = MMF_GRS;
    
end    

%Average Results over Channels
MMF_RS_average = mean(MMF_RS_set);
MMF_NoRS_average = mean(MMF_NoRS_set);
MMF_RS_cluster_average = mean(MMF_RS_cluster_set);
MMF_HRS_average = mean(MMF_HRS_set);
MMF_GRS_average = mean(MMF_GRS_set);

%Plot results
plot(PACs,MMF_RS_average,'-sr','Linewidth',1.2)
hold on
grid on
plot(PACs,MMF_NoRS_average,'-sm','Linewidth',1.2)
plot(PACs,MMF_RS_cluster_average,'-sg','Linewidth',1.2)
plot(PACs,MMF_HRS_average,'-sb','Linewidth',1.2)
plot(PACs,MMF_GRS_average,'-sk','Linewidth',1.2)

xlabel('Per antenna feed power (W)')
ylabel('Rate (bit/s/Hz)')
title('Max-Min Fair Multigroup Multicast')
legend('1-layer RS, perfect','NoRS, perfect','RS per cluster, perfect','HRS, perfect', 'GRS, perfect',...
'Location','northwest')
