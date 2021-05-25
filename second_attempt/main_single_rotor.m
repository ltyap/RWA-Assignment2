%% single rotor calculations
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

%% other parameters
a_wake_bem = 0.2602;
N = 20; % nr of blade segments
spacing = 1; % 1- uniform, 0-cosine
L = 0;    % distance between 2 rotors
sec_rot=0;  % is there a second rotor
phase_diff = 0;
Nrotations = 10; %number of rotations done by the rotor (length of wake)

[CT,CP,CQ, results] = lifting_line_loop(rotor,windvel,N,spacing,L, sec_rot, a_wake_bem, Nrotations);
filename = "results/results_llt_N"+N+".mat";
save(filename,'results');

plot_BEM_vs_LLT(windvel, rotor.TSR, rotor.Radius, rotor.NBlades, rho);

%% influence of different parameters

% assumed convection speed of wake
a_wake = [0,0.25,a_wake_bem, 0.5,0.75];
filenames={};
for i=1:length(a_wake)
    [CT_a(i),CP_a(i),CQ_a(i), results] = lifting_line_loop(rotor, windvel, N,spacing,L, sec_rot, a_wake(i), Nrotations, phase_diff)
    filename = "results/results_llt_N"+N+"a_"+a_wake(i)+".mat";
    filenames{i} = filename;
    save(filename,'results');
end
var_change = "a_{wake}";% variable that changes - for plot legend
plot_compare(filenames, a_wake, var_change, N, rotor, rho, windvel);

% discretization of the blade (constant, cosine) azimuthal discretization 
% (number of wake segments per rotation)
spacing=[0,1];
filenames={};
for i=1:length(spacing)
    [CT_cos(i),CP_cos(i),CQ_cos(i),results] = lifting_line_loop(rotor, windvel, N,spacing,L, sec_rot, a_wake_bem, Nrotations, phase_diff);
    filename = "results/results_llt_N"+N+"cos_"+spacing(i)+".mat";
    filenames{i} = filename;
    save(filename,'results');
end
var_change = "spacing";% variable that changes - for plot legend
plot_compare(filenames, spacing, var_change, N, rotor, rho, windvel);


% length of the wake (number of rotations), including convergence of the solution with wake length
spacing=1;
Nrotations = [0.1, 1, 5,10,15];
for i=1:length(Nrotations)
    [CT_rot(i),CP_rot(i),CQ_rot(i),results, conv(i)] = lifting_line_loop(rotor, windvel, N,spacing,L, sec_rot, a_wake_bem, Nrotations(i), phase_diff);
    filename = "results/results_llt_N"+N+"rot_"+Nrotations(i)+".mat";
    filenames{i} = filename;
    save(filename,'results');
end
var_change = "N_{rotations}";% variable that changes - for plot legend
plot_compare(filenames, Nrotations, var_change, N, rotor, rho, windvel);
