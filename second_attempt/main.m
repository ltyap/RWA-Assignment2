%% Assignment 2 of RWA - lifting line model
clc
clear all
close all
%% define the blade geometry - same as in BEM model
% blade discretisation
TipLocation_R =  1;%non-dimensional]
RootLocation_R =  0.2;%[non-dimensional]
N = 8; % number of segments

% spacing: 1 for uniform, 0 for cosine
spacing = 1;
[dist] = RadialSpacing(N, TipLocation_R, RootLocation_R, spacing);
r_R = dist.r_R;
%% flow conditions - same as in BEM model
windvel = [10,0,0];      % freestream velocity [m/s]
altitude = 0;%[km]
[~, ~, pinf, rho] = atmosisa(altitude);
%% rotor parameters
TSR = 8;        % tip speed ratios we want to calculate
Radius = 50;    % blade radius(length) [m]
NBlades = 3;    % number of blades
a_wake = 0.2602;   % ,should be average induction at the rotor, from BEM
Omega = windvel(1)*TSR/Radius;
Nrotations = 0.5;%for the wake
% theta_array = linspace(0,Nrotations*2*pi);%Omega*t, where t is the time
theta_array = [0:pi/10:2*pi*Nrotations];%Omega*t, where t is the time
% Lw_D:  wake length in diameters downstream
Lw_D = max(theta_array)/Omega*windvel(1)*(1-a_wake)/(2*Radius);% [-]
%% Main 
RotorWakeSystem = vortex_system(r_R, Radius, TSR/(1-a_wake), theta_array, NBlades);

[InfluenceMatrix] = InfluenceMatrix(RotorWakeSystem, NBlades);
[a, aline, r_R_cp, Fnorm, Ftan, Gamma_temp]= solveSystem(InfluenceMatrix, RotorWakeSystem, Radius, Omega, windvel);
[CT, CP, CQ] = CT_CPcalculations(Fnorm, Ftan, windvel(1), r_R, r_R_cp, Omega, Radius, NBlades);

% %%%%%% SHOULD BE %%%%%%%%%
% % CT = 0.6553;
% % CQ = 0.2801;
% % CP = 0.4482;