%% all dual rotor calculations
clc
clear all
close all
%% define the blade geometry - same as in BEM model
% blade discretisation
rotor.TipLocation_R =  1;     % non-dimensional
rotor.RootLocation_R =  0.2;  % non-dimensional

%% rotor parameters
rotor.TSR = 8;        % tip speed ratios we want to calculate
rotor.Radius = 50;    % blade radius(length) [m]
rotor.NBlades = 3;    % number of blades

%% flow conditions - same as in BEM model
windvel = [10,0,0]; % freestream velocity [m/s]
altitude = 0;   %[km]
[~, ~, pinf, rho] = atmosisa(altitude);

%% operational conditions
D = 2*rotor.Radius;
N = 20;
Nwake = 201;    % nr of wake discretization
spacing=1; % 1- uniform, 0-cosine 
sec_rot=1;  % is there a second rotor
a_wake_bem = 0.2602;
Nrotations = 1; %number of rotations done by the rotor (length of wake)
L = [1,2,5,1000]*D; % distance between rotors
% L=[1,2]*D;
phase_diff = 0; % in degrees

%% calculations
names={};
for i=1:length(L)
    [CT(i,:),CP(i,:),CQ(i,:), results] = lifting_line_loop(rotor,windvel,N,spacing,L(i),sec_rot, a_wake_bem, Nrotations, Nwake, phase_diff);
    filename = "results/results_dual_llt_N"+N+"_phase_"+phase_diff+"_L"+L(i)+".mat";
    filenames{i} = filename;
    save(filename,'results');
end
var_change = "L";% variable that changes - for plot legend
plot_compare(filenames, L, var_change, N, rotor, rho, windvel);
