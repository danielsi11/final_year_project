clear

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

%Satellite System Parameter
theta_3dB = 0.4*pi/180;
height = 35786*10^3;

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

beam = horzcat(beam_temp+shift1,beam_temp+shift2,beam_temp+shift3);

%Plot beam centres
p(1) = plot(beam(1,:),beam(2,:),'r*');
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
p(2) = plot(user(1,:),user(2,:),'b*');

%Change group clustering
beam = beam(:,[1 4 7 2 5 8 3 6 9]);

%Plot beam boundaries
for i_beam = 1:3
    x = beam(1,i_beam);
    y = beam(2,i_beam);
    th = 0:0.01:2*pi;
    xunit = beam_radius*cos(th) + x;
    yunit = beam_radius*sin(th) + y;
    p(3) = plot(xunit, yunit,'m','Linewidth',2);
end

for i_beam = 4:6
    x = beam(1,i_beam);
    y = beam(2,i_beam);
    th = 0:0.01:2*pi;
    xunit = beam_radius*cos(th) + x;
    yunit = beam_radius*sin(th) + y;
    p(4) = plot(xunit, yunit,'g','Linewidth',2);
end

for i_beam = 7:9
    x = beam(1,i_beam);
    y = beam(2,i_beam);
    th = 0:0.01:2*pi;
    xunit = beam_radius*cos(th) + x;
    yunit = beam_radius*sin(th) + y;
    p(5) = plot(xunit, yunit,'y','Linewidth',2);
end

hold off
title('Beam Pattern and User Terminal Location')
xlim([-900000 900000])
ylim([-900000 900000])
legend(p,'Beam Center','User Terminal','Group Cluster 1','Group Cluster 2','Group Cluster 3','Location','northwest')
