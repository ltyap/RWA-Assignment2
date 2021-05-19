%% TO DO
% add control points from all blades, not just one

% Assignment 2 of RWA - lifting line model
clc
clear all
close all

[polar_alpha, polar_cl, polar_cd]=import_polars();
%% define the blade geometry - same as in BEM model
% blade discretisation
TipLocation_R =  1;
RootLocation_R =  0.2;
N = 10; % number of segments

% uniform distribution
% delta_r_R = 0.01;
% r_R = [0.2:delta_r_R:1];

[r_R, theta, b] = RadialSpacing(N, TipLocation_R, RootLocation_R);% cosine distribution for now
[chord_distribution, twist_distribution] = BladeGeometry(r_R);
%% flow conditions - same as in BEM model
Uinf = 10;      % freestream velocity [m/s]
rho = 1.225;    % air density
pinf = 101300;  % atmospheric pressure
TSR = 8;        % tip speed ratios we want to calculate
Radius = 50;    % blade Radius [m]
NBlades = 3;    % number of blades
a_wake = 0.2602;   % ,should be average induction at the rotor, from BEM
Uwake = Uinf*(1-a_wake);

Omega = Uinf*TSR/Radius;
Lw_D = 0.5;    % wake length in diameters downstream

[cp, bound, trail]=vortex_system(Omega, NBlades, N, r_R, Radius, twist_distribution, chord_distribution, Lw_D, Uwake);

%% determine induced velocity per unit strength of circulation (for Gamma=1)
for i=1:length(cp.x) % take all control points
    local_cp = [cp.x(i), cp.y(i), cp.z(i)];
    for n=1:NBlades
        for j=1:N
            induced_vel=[0,0,0];
            % bound vortex influence
            point1 = [bound.x((n-1)*N+j), bound.y((n-1)*N+j), bound.z((n-1)*N+j)]; % coordinates of start of bound vortex
            point2 = [bound.x((n-1)*N+j+1), bound.y((n-1)*N+j+1), bound.z((n-1)*N+j+1)]; % coordinates of end of bound vortex
            temp = induced_v_from_vortex(1, point1, point2, local_cp);
            induced_vel=induced_vel+temp;

            % first trailing vortex
            temp = induced_v_from_vortex(1,[trail.x((n-1)*N+j,1),trail.y((n-1)*N+j,1),trail.z((n-1)*N+j,1)] ,point1, local_cp);
            induced_vel = induced_vel+temp;

            temp = induced_v_from_vortex(1, point2, [trail.x((n-1)*N+j+1,1),trail.y((n-1)*N+j+1,1),trail.z((n-1)*N+j+1,1)], local_cp);
            induced_vel = induced_vel+temp;  

            for k=2:length(trail.x((n-1)*N+j,:))-1 % for the whole length of wake - positive circulation
                temp = induced_v_from_vortex(1,[trail.x((n-1)*N+j,k+1),trail.y((n-1)*N+j,k+1),trail.z((n-1)*N+j,k+1)] ,[trail.x((n-1)*N+j,k),trail.y((n-1)*N+j,k),trail.z((n-1)*N+j,k)], local_cp);
                induced_vel = induced_vel+temp;
            end

            for k=2:length(trail.x((n-1)*N+j,:))-1 % for the whole length of wake - negative circulation
                temp = induced_v_from_vortex(1,[trail.x((n-1)*N+j+1,k),trail.y((n-1)*N+j+1,k),trail.z((n-1)*N+j+1,k)] ,[trail.x((n-1)*N+j+1,k+1),trail.y((n-1)*N+j+1,k+1),trail.z((n-1)*N+j+1,k+1)], local_cp);
                induced_vel = induced_vel+temp;
            end

            % write into influence matrix - for unit gamma
            Influence_u(i,(n-1)*N+j) = induced_vel(1);
            Influence_v(i,(n-1)*N+j) = induced_vel(2);
            Influence_w(i,(n-1)*N+j) = induced_vel(3);
        end
    end
end
 
%% Calculate
%why do we also have r_R as an output?
[a, aline, r_R, Fnorm, Ftan, Gamma] = solveGamma(Uinf, N, Radius, cp, Influence_u, Influence_v, Influence_w, Omega, polar_alpha, polar_cl, polar_cd, chord_distribution);

%% Post processing
% [r_R, theta, b] = RadialSpacing(N, TipLocation_R, RootLocation_R);
temp = RootLocation_R + flip(b*(cos(theta)+1));  
dr = temp(2:end)-temp(1:end-1); % lengths of sections
dr = repmat(dr,1,3);

CT_sections = Fnorm*NBlades.*dr/(0.5*Uinf^2*pi*Radius^2); % local thrust coefficient
CQ_sections = Ftan.*r_R*NBlades.*dr*Radius/(0.5*Uinf^3*pi*Radius^2);    % local torque coefficient
CP_sections = Ftan.*r_R*NBlades.*dr*Radius*Omega/(0.5*Uinf^3*pi*Radius^2); % local power coefficient  P=omega*Q
Gamma = Gamma/(pi*Uinf^2/(NBlades*Omega)); % non-dimensionalize circulation

CT = sum(dr.*Fnorm*NBlades/(0.5*Uinf^2*pi*Radius^2));     % total thrust coefficient
CP = sum(dr'.*Ftan.*r_R*NBlades*Radius*Omega/(0.5*Uinf^3*pi*Radius^2));    % total power coefficient

CT_all = sum(CT_sections, 'all');  % write to array of thrust coefficients
CQ_all = sum(CQ_sections, 'all');  % write to array of torque coefficients
CP_all = sum(CP_sections, 'all');  % write to array of power coefficients

% %%%%%% SHOULD BE %%%%%%%%%
% % CT = 0.6553;
% % CQ = 0.2801;
% % CP = 0.4482;

