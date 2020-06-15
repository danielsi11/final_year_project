clc; clear all; clf

%Number of Transmitter Antennas
N = 7; 

%Number of Groups
M = 7;

%Number of Users per Group
rho = 2;

%Total Number of Users
K = rho*M;

%Number of Users in Each Group
Gm = repmat(rho,1,M);

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
theta = (0:N-2)*2*pi/(N-1);
beam_spacing = ((beam_radius*3)/2)/cos(pi/(N-1));

beam = [beam_spacing*cos(theta); beam_spacing*sin(theta); zeros(1,N-1)];
beam = horzcat([0;0;0], beam);

%Range of CSIT Error Scaling Factor
alphas = 0.6:0.2:0.8;

%Channel Sample Size
S = 1000;

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

            P = PAC*N; %Total Transmit Power
            P_e = P^(-alpha); %Error Variance          
            
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
                
                %Channel Error
                H_error(:,i_user) = sqrt(P_e)/sqrt(2)*(randn(N,1)+1i*randn(N,1));

            end

            %Channel Estimate
            H_est = H - H_error;

            if i_alpha == length(alphas)+1 %Perfect CSIT
                MMF_RS(i_PAC) = RS_perfect_rate(Gm,H,PAC,tolerance);
                MMF_NoRS(i_PAC) = NoRS_perfect_rate(Gm,H,PAC,tolerance); 

            else %Imperfect CSIT
                MMF_RS(i_PAC) = RS_imperfect_rate(Gm,H_est,PAC,tolerance,S,P_e);
                MMF_NoRS(i_PAC) = NoRS_imperfect_rate(Gm,H_est,PAC,tolerance,S,P_e); 
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

plot(PACs,MMF_RS_average(1,:,1),'--sr','Linewidth',1.2)
plot(PACs,MMF_NoRS_average(1,:,1),'--sb','Linewidth',1.2)

xlabel('Per antenna feed power (W)')
ylabel('Rate (bit/s/Hz)')
title('Max-Min Fair Multigroup Multicast')
legend('RS, perfect','NoRS, perfect','RS, imperfect 0.8','NoRS, imperfect 0.8',...
'RS, imperfect 0.6','NoRS, imperfect 0.6','Location','northwest')
