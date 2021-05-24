%% single rotor calculations
clc
clear all
close all

a_wake_bem = 0.2602;
N = 20;
spacing = 1; % 1- uniform, 0-cosine
L=0;    % distance between 2 rotors
sec_rot=0;  % is there a second rotor

% [CT,CP,CQ] = lifting_line_loop(N,spacing,L, sec_rot,a_wake_bem);

% plot_BEM_vs_LLT;

%% influence of different parameters

% assumed convection speed of wake
a_wake = [0,0.25,a_wake_bem, 0.5,0.75];
for i=1:length(a_wake)
    [CT_a(i),CP_a(i),CQ_a(i), results_a] = lifting_line_loop(N,spacing,L, sec_rot,a_wake(i));
    filename = "results_llt_N"+N+"a_"+a_wake(i)+".mat";;
    save(filename,'results_a');
end
plot_compare(a_wake,N);

% discretization of the blade (constant, cosine) azimuthal discretization 
% (number of wake segments per rotation)
spacing=[0,1];
for i=1:length(spacing)
    [CT_cos(i),CP_cos(i),CQ_cos(i)] = lifting_line_loop(N,spacing(i),L, sec_rot,a_wake_bem);
end
plot_compare;
% length of the wake (number of rotations), including convergence of the solution with wake length
spacing=1;
Nrotations = [0.1, 1, 5,10,15];
for i=1:length(L_D)
    [CT_rot(i),CP_rot(i),CQ_rot(i),results_rot, conv(i)] = lifting_line_loop(N,spacing,L, sec_rot,a_wake_bem, Nrotations);
end
