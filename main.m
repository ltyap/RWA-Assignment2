%% TO DO
% add control points from all blades, not just one

% Assignment 2 of RWA - lifting line model
clc
clear all
close all

[polar_alpha, polar_cl, polar_cd]=import_polars();

%% define the blade geometry - same as in BEM model
% uniform distribution
% delta_r_R = 0.01;
% r_R = [0.2:delta_r_R:1];

% blade shape
TipLocation_R =  1;
RootLocation_R =  0.2;

% cosine distribution 
N = 31; % number of segments
theta = linspace(0,pi,N+1);
b = (TipLocation_R-RootLocation_R)/2;
r_R = RootLocation_R + flip(b*(cos(theta)+1)); % cosine spacing

pitch = 2; % degrees
chord_distribution = 3*(1-r_R)+1; % meters
twist_distribution = -14*(1-r_R)+pitch; % degrees

%% flow conditions - same as in BEM model
Uinf = 10;      % freestream velocity
rho = 1.225;    % air density
pinf = 101300;  % atmospheric pressure
TSR = 8;        % tip speed ratios we want to calculate
Radius = 50;    % blade Radius
NBlades = 3;    % number of blades
a_wake = 0.2602;   % ,should be average induction at the rotor, from BEM
Uwake = Uinf*(1-a_wake);

Omega = Uinf*TSR/Radius;
Lw_D = 0.2;    % wake length in diameters downstream

theta_array = [0:2*pi/NBlades:2*pi*(1-1/NBlades)]; % azimuthal positions of blades, assuming the first one is at 0

[cp, bound, trail]=vortex_system(Omega, NBlades, theta_array,N, r_R, Radius, twist_distribution, chord_distribution, Lw_D, Uwake)

%% determine induced velocity per unit strength of circulation (for Gamma=1)
% for i=1:length(cp(:,1))
%     for j=1:length(r_R)-1
%         induced_vel=[0,0,0];
%         
%         point1 = [0,chord_distribution(j)/4,r_R(j)*Radius]; % coordinates of start of bound vortex
%         point2 = [0,chord_distribution(j+1),r_R(j+1)*Radius]; % coordinates of end of bound vortex
%         temp = induced_v_from_vortex(1,point1, point2, cp(i,:));
%         induced_vel=induced_vel+temp;
%         
%         % first trailing vortex
%         temp = induced_v_from_vortex(1,[trail_points.x(j,1),trail_points.y(j,1),trail_points.z(j,1)] ,point1, cp(i,:));
%         induced_vel = induced_vel+temp;
%         
%         temp = induced_v_from_vortex(1, point2, [trail_points.x(j+1,1),trail_points.y(j+1,1),trail_points.z(j+1,1)], cp(i,:));
%         induced_vel = induced_vel+temp;        
%         
%         for k=2:length(trail_points.x(i,:))-1 % for the whole length of wake - positive circulation
%             temp = induced_v_from_vortex(1,[trail_points.x(j,k+1),trail_points.y(j,k+1),trail_points.z(j,k+1)] ,[trail_points.x(j,k),trail_points.y(j,k),trail_points.z(j,k)], cp(i,:));
%             induced_vel = induced_vel+temp;
%         end
%         
%         for k=2:length(trail_points.x(i,:))-1 % for the whole length of wake - negative circulation
%             temp = induced_v_from_vortex(1,[trail_points.x(j+1,k),trail_points.y(j+1,k),trail_points.z(j+1,k)] ,[trail_points.x(j+1,k+1),trail_points.y(j+1,k+1),trail_points.z(j+1,k+1)], cp(i,:));
%             induced_vel = induced_vel+temp;
%         end
%         
%         % write into influence matrix - for unit gamma
%         Influence_u(i,j) = induced_vel(1);
%         Influence_v(i,j) = induced_vel(2);
%         Influence_w(i,j) = induced_vel(3);
%     end
% end
% 
% %% check if the system if set up correctly
% figure('visible', 'on')
% hold on
% plot3( zeros(1,2*(N+1)), [zeros(1,N+1),flip(chord_distribution)], [r_R*Radius,flip(r_R*Radius)])% blade outline
% scatter3(cp(:,1), cp(:,2), cp(:,3), 15, 'r', 'filled')  % control points
% plot3(zeros(1,N+1),chord_distribution/4,r_R*Radius,'-+', 'Color', 'k')   % bound vortices
% for i=1:N+1
%     plot3([0, trail_points.x(i,:)], [chord_distribution(i)/4,trail_points.y(i,:)],[r_R(i)*Radius,trail_points.z(i,:)], 'k')
% end
% % axis equal
% xlabel("x")
% ylabel("y")
% zlabel("z")
% hold off
% 
% %% Calculate
% %why do we also have r_R as an output?
% [a, aline, r_R, Fnorm, Ftan, Gamma] = solveGamma(Uinf, N, Radius, cp, Influence_u, Influence_v, Influence_w, Omega, polar_alpha, polar_cl, polar_cd, cp_twist);
% 
% %% Post processing
% temp = RootLocation_R + flip(b*(cos(theta)+1));  
% dr = temp(2:end)-temp(1:end-1); % lengths of sections
% 
% CT_sections = Fnorm*NBlades.*dr/(0.5*Uinf^2*pi*Radius^2); % local thrust coefficient
% CQ_sections = Ftan.*r_R*NBlades.*dr*Radius/(0.5*Uinf^3*pi*Radius^2);    % local torque coefficient
% CP_sections = Ftan.*r_R*NBlades.*dr*Radius*Omega/(0.5*Uinf^3*pi*Radius^2); % local power coefficient  P=omega*Q
% Gamma = Gamma/(pi*Uinf^2/(NBlades*Omega)); % non-dimensionalize circulation
% 
% CT = sum(dr.*Fnorm*NBlades/(0.5*Uinf^2*pi*Radius^2));     % total thrust coefficient
% CP = sum(dr'.*Ftan.*r_R*NBlades*Radius*Omega/(0.5*Uinf^3*pi*Radius^2));    % total power coefficient
% 
% CT_all = sum(CT_sections, 'all');  % write to array of thrust coefficients
% CQ_all = sum(CQ_sections, 'all');  % write to array of torque coefficients
% CP_all = sum(CP_sections, 'all');  % write to array of power coefficients
% 
% %%%%%% SHOULD BE %%%%%%%%%
% % CT = 0.6553;
% % CQ = 0.2801;
% % CP = 0.4482;

