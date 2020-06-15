clc; clear all; clf;

%Number of Transmitter Antennas
N = 9; 

%Number of Gateways
L = 3;

%Number of Transmitter Antennas per Feed Cluster
B = N/L; 

%Number of Groups
M = 9;

%Number of Users per Group
rho = 2;

%Total Number of Users
K = rho*M;

%Number of Users in Each Group
Gm = repmat(rho,B,L);

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

beam(:,:,1) = beam_temp + shift1;
beam(:,:,2) = beam_temp + shift2;
beam(:,:,3) = beam_temp + shift3;

%Range of CSIT Error Scaling Factor
alphas = 0.6:0.2:0.8;

%Channel Sample Size
S = 1000;

%Feeder Link Interference Level
delta = 0.8;

%Simulate MMF rate
for i_channel = 1:N_channel

    for i_alpha = 1:length(alphas)+1

        if i_alpha ~= length(alphas)+1
            alpha = alphas(i_alpha);
        end

        for i_PAC = 1:length(PACs)   

            rng('default');
            rng(i_channel+1);

            PAC = PACs(i_PAC);

            P = PAC*B*L; %Total Transmit Power
            P_e = P^(-alpha); %Error Variance

            offset_group = 0;
            offset_cluster = 0;
            i_group = 1;
            i_cluster = 1;
            for i_user = 1:K

                %User Distribution
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

                theta_rng = 2*pi*rand;
                user(:,i_user) = beam(:,i_group,i_cluster) + beam_radius*rand*[cos(theta_rng); sin(theta_rng); 0];

            end
            
            for i_user = 1:K

                %User Link Channel Model
                X_k_dB(:,i_user) = exp(rain_mu + rain_sigma*randn*ones(B,1));
                X_k(:,i_user) = 10.^(X_k_dB(:,i_user)/20);

                phi_k(:,i_user) = 2*pi*rand*ones(B,1);

                q(:,i_user) = X_k(:,i_user).^(-0.5).*exp(-1i*phi_k(:,i_user));

                d_k(i_user) = norm(satellite - user(:,i_user));

                user_sat = user(:,i_user) - satellite;

                for i_cluster = 1:L
                    for i_beam = 1:B
                        beam_sat = beam(:,i_beam,i_cluster) - satellite;
                        theta_nk(i_beam,i_user,i_cluster) = acos((beam_sat'*user_sat)/(norm(beam_sat)*norm(user_sat)));
                    end

                    u_nk(:,i_user,i_cluster) = 2.07123*sin(theta_nk(:,i_user,i_cluster))/sin(theta_3dB);

                    G_nk(:,i_user,i_cluster) = G_max*(besselj(1,u_nk(:,i_user,i_cluster))./(2*u_nk(:,i_user,i_cluster)) + ...
                    36*besselj(3,u_nk(:,i_user,i_cluster))./(u_nk(:,i_user,i_cluster).^3)).^2;

                    b(:,i_user,i_cluster) = sqrt(G_R*G_nk(:,i_user,i_cluster))./(4*pi*(d_k(i_user)/wavelength)*sqrt(Boltz*T_sys*B_w));

                    H(:,i_user,i_cluster) = b(:,i_user,i_cluster).*q(:,i_user);

                    %Channel Error
                    H_error(:,i_user,i_cluster) = sqrt(P_e)/sqrt(2)*(randn(B,1)+1i*randn(B,1));

                end     
            end 

            %Feeder Link Channel Model
            for i_gateway = 1:L
                X_l_dB(i_gateway) = exp(rain_mu + rain_sigma*randn);
                X_l(i_gateway) = 10^(X_l_dB(i_gateway)/20);

                phi_l(i_gateway) = 2*pi*rand;

                q_l(i_gateway) = X_l(i_gateway)^(-0.5)*exp(-1i*phi_l(i_gateway));

                for i_cluster = 1:L
                    if i_cluster == i_gateway
                        F(:,:,i_cluster,i_gateway) = eye(B);
                    else
                        F(:,:,i_cluster,i_gateway) = delta*ones(B,B);
                    end

                    G(:,:,i_cluster,i_gateway) = q_l(i_gateway)*F(:,:,i_cluster,i_gateway);

                end

            end

            if i_alpha == length(alphas)+1 %Perfect CSIT
                MMF_RS(i_PAC) = RS_perfect_rate(Gm,H,PAC,tolerance,G);
                MMF_NoRS(i_PAC) = NoRS_perfect_rate(Gm,H,PAC,tolerance,G);               

            else %Imperfect CSIT
                MMF_RS(i_PAC) = RS_imperfect_rate(Gm,H,H_error,PAC,tolerance,S,P_e,G);
                MMF_NoRS(i_PAC) = NoRS_imperfect_rate(Gm,H,H_error,PAC,tolerance,S,P_e,G); 
            end 

        end

        MMF_RS_set(i_channel,:,i_alpha) = MMF_RS;
        MMF_NoRS_set(i_channel,:,i_alpha) = MMF_NoRS;

    end
    
end

%Average Results over Channels
MMF_RS_average = mean(MMF_RS_set);
MMF_NoRS_average = mean(MMF_NoRS_set);

%Plot results
plot(PACs,MMF_RS_average(1,:,3),'-sr','Linewidth',1.2)
hold on
grid on
plot(PACs,MMF_NoRS_average(1,:,3),'-sb','Linewidth',1.2)

plot(PACs,MMF_RS_average(1,:,2),'--or','Linewidth',1.2)
plot(PACs,MMF_NoRS_average(1,:,2),'--ob','Linewidth',1.2)

plot(PACs,MMF_RS_average(1,:,1),'-.vr','Linewidth',1.2)
plot(PACs,MMF_NoRS_average(1,:,1),'-.vb','Linewidth',1.2)

xlabel('Per antenna feed power (W)')
ylabel('Rate (bit/s/Hz)')
title('Max-Min Fair Multigroup Multicast')
legend('RS, perfect','NoRS, perfect','RS, imperfect 0.8','NoRS, imperfect 0.8',...
'RS, imperfect 0.6','NoRS, imperfect 0.6','Location','northwest')
