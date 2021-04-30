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
TSR_array =[8];    %TSR_array =[6,8,10]; % tip speed ratios we want to calculate
Radius = 50;    % blade Radius
NBlades = 3;    % number of blades

TSR = TSR_array(1);
Omega = Uinf*TSR/Radius;
results =zeros(length(r_R)-1,9);

% determine control point locations - at middle of segments
c = (chord_distribution(1:end-1)+chord_distribution(2:end))/2;  % local chord
cp(:,1) = c/4; % x-coordinate
cp(:,2) = (r_R(1:end-1)+r_R(2:end))/2*Radius;     % y-coordinate
cp(:,3) = zeros(N,1);     % z-coordinate

% twist at control points
cp_twist = interp1(r_R, twist_distribution, cp(:,2)/Radius);

% determine induced velocity per unit strength of circulation (for Gamma=1)
% influence matrices
trail_vortex = 1e6; % length of trailing edge vortex, for now we assume straight, infinity approximated by 1e6
for i=1:length(cp(:,1))
    % we know this will be zero for this case (rectangular wing), but not
    % for others
    for j=1:length(r_R)-1
        induced_vel=[0,0,0];
        point1 = [chord_distribution(j)/4,r_R(j)*Radius,0]; % coordinates of start of bound vortex
        point2 = [chord_distribution(j+1),r_R(j+1)*Radius,0]; % coordinates of end of bound vortex
        temp = induced_v_from_vortex(1,point1, point2, cp(i,:));
        induced_vel=induced_vel+temp;
        
        % add influence of two trailing vortices
        temp = induced_v_from_vortex(1,point1+[trail_vortex,0,0], point1, cp(i,:));
        induced_vel=induced_vel+temp;
        
        temp = induced_v_from_vortex(1,point2, point2+[trail_vortex,0,0], cp(i,:));
        induced_vel=induced_vel+temp;
        
        % write into influence matrix - for unit gamma
        Influence_u(i,j) = induced_vel(1);
        Influence_v(i,j) = induced_vel(2);
        Influence_w(i,j) = induced_vel(3);
    end
end

% why do we also have r_R as an output?
[a, aline, r_R, Fnorm, Ftan, Gamma] = solveGamma(Uinf, N, Radius, cp, Influence_u, Influence_v, Influence_w, Omega, polar_alpha, polar_cl, polar_cd, cp_twist);

%% Post processing
temp = RootLocation_R + flip(b*(cos(theta)+1));  
dr = temp(2:end)-temp(1:end-1); % lengths of sections

CT_sections = Fnorm*NBlades.*dr'/(0.5*Uinf^2*pi*Radius^2); % local thrust coefficient
CQ_sections = Ftan.*r_R*NBlades.*dr'*Radius/(0.5*Uinf^3*pi*Radius^2);    % local torque coefficient
CP_sections = Ftan.*r_R*NBlades.*dr'*Radius*Omega/(0.5*Uinf^3*pi*Radius^2); % local power coefficient  P=omega*Q
Gamma = Gamma/(pi*Uinf^2/(NBlades*Omega)); % non-dimensionalize circulation

CT = sum(dr.*Fnorm*NBlades/(0.5*Uinf^2*pi*Radius^2));     % total thrust coefficient
CP = sum(dr'.*Ftan.*r_R*NBlades*Radius*Omega/(0.5*Uinf^3*pi*Radius^2));    % total power coefficient

CT_all = sum(CT_sections, 'all');  % write to array of thrust coefficients
CQ_all = sum(CQ_sections, 'all');  % write to array of torque coefficients
CP_all = sum(CP_sections, 'all');  % write to array of power coefficients

