clear

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

%Satellite System Parameter
theta_3dB = 0.4*pi/180;
height = 35786*10^3;

%Beam Location
beam_radius = height*tan(theta_3dB); 
theta = (0:N-2)*2*pi/(N-1);
beam_spacing = ((beam_radius*3)/2)/cos(pi/(N-1));

beam = [beam_spacing*cos(theta); beam_spacing*sin(theta); zeros(1,N-1)];

beam = horzcat([0;0;0], beam);

%Plot beam centres
plot(beam(1,:),beam(2,:),'r*');
hold on

%User Distribution
offset = 0;
i_group = 1;
for i_user = 1:K
    if(i_user > Gm(i_group) + offset)
        offset = offset + Gm(i_group);
        i_group = i_group + 1;
    end
    theta_rng = 2*pi*rand;
    user(:,i_user) = beam(:,i_group) + beam_radius*rand*[cos(theta_rng); sin(theta_rng); 0];
    
end    

%Plot user terminals
plot(user(1,:),user(2,:),'b*');

%Plot beam boundaries
for i_beam = 1:N
    x = beam(1,i_beam);
    y = beam(2,i_beam);
    th = 0:0.01:2*pi;
    xunit = beam_radius*cos(th) + x;
    yunit = beam_radius*sin(th) + y;
    plot(xunit, yunit,'m','Linewidth',2);
end
hold off
title('Beam Pattern and User Terminal Location')
legend('Beam Center','User Terminal','Beam Boundary','Location','northwest')
