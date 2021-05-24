function [CT,CP,CQ] = lifting_line_loop(N,spacing,L, sec_rot)
%% Assignment 2 of RWA - lifting line model
%% define the blade geometry - same as in BEM model
% blade discretisation
TipLocation_R =  1;     % non-dimensional
RootLocation_R =  0.2;  % non-dimensional
%N = 20;  % number of segments

% spacing: 1 for uniform, 0 for cosine
% spacing = 1;
[r_R] = RadialSpacing(N, TipLocation_R, RootLocation_R, spacing);
%i think it works for: N = 11, error tolerance = 1e-4 and Nrotations = 15,
%Separation distance  = 4D
%% flow conditions - same as in BEM model
windvel = [10,0,0];% freestream velocity [m/s]
altitude = 0;%[km]
[~, ~, pinf, rho] = atmosisa(altitude);
%% rotor parameters
TSR = 8;        % tip speed ratios we want to calculate
Radius = 50;    % blade radius(length) [m]
NBlades = 3;    % number of blades
a_wake = 0.2602;   % average induction at the rotor, from BEM
Omega = norm(windvel)*TSR/Radius;

Nrotations = 10;   % for the wake
theta_array = [0:pi/10:2*pi*Nrotations];%Omega*t, where t is the time
% Lw_D:  wake length in diameters downstream
Lw_D = max(theta_array)/Omega*norm(windvel)*(1-a_wake)/(2*Radius); % [-]

%% LLT calculations
RotorWakeSystem = vortex_system(r_R, Radius, TSR/(1-a_wake), theta_array, NBlades,sec_rot,L);

[InfluenceMatrix] = influence_matrix(RotorWakeSystem);
[a, aline, r_R_cp, Fnorm, Ftan, GammaNew, alpha, inflow]= solveSystem(InfluenceMatrix, RotorWakeSystem, Radius, Omega, windvel);
%QUICKFIX
if sec_rot==1
    r_R_cp = [r_R_cp(1:N*NBlades);r_R_cp(1:N*NBlades)];
end
[CT, CP, CQ, ct, cp, cq] = CT_CPcalculations(N,Fnorm, Ftan, windvel(1), r_R, r_R_cp, Omega, Radius, NBlades, sec_rot);

% %%%%%% FROM BEM %%%%%%%%%
% % CT = 0.6553;
% % CQ = 0.2801;
% % CP = 0.4482;

% load('tsr6.mat');
% load('tsr8.mat');
% load('tsr10.mat');

plotting_func(sec_rot,windvel,Radius, N, NBlades, Omega, a, aline, r_R_cp, ct,cp,cq, GammaNew, alpha, inflow);

%% saving data to file
results_llt = [r_R_cp,a,aline,ct,cp,cq,GammaNew,alpha',inflow', Fnorm, Ftan];
if sec_rot==1
    filename = "results_llt_dual_N"+N+"L_"+L+".mat";
else
    filename = "results_llt_N"+N+".mat";
end
save(filename,'results_llt');
end