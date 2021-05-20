%% Assignment 2 of RWA - lifting line model
clc
clear all
close all

% [polar_alpha, polar_cl, polar_cd]=import_polars();
%% define the blade geometry - same as in BEM model
% blade discretisation
TipLocation_R =  1;%non-dimensional]
RootLocation_R =  0.2;%[non-dimensional]
N = 10; % number of segments

% cosine distribution for now
[r_R, ~, ~] = RadialSpacing(N, TipLocation_R, RootLocation_R);
% [chord_distribution, twist_distribution] = BladeGeometry(r_R);
%% flow conditions - same as in BEM model
Uinf = [10,0,0];      % freestream velocity [m/s]
altitude = 0;%[km]
[~, ~, pinf, rho] = atmosisa(altitude);

% rotor parameters
TSR = 8;        % tip speed ratios we want to calculate
Radius = 50;    % blade radius(length) [m]
NBlades = 3;    % number of blades
a_wake = 0.2602;   % ,should be average induction at the rotor, from BEM
% Uwake = Uinf(1)*(1-a_wake);
Omega = Uinf(1)*TSR/Radius;
Nrotations = 1;%for the wake
theta_array = linspace(0,Nrotations*2*pi);%Omega*t, where t is the time
RotorWakeSystem = vortex_system(r_R, Radius, TSR/(1-a_wake), theta_array, NBlades);

InfluenceMatrix();


