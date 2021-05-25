function [CT,CP,CQ, results,conv] = lifting_line_loop(rotor, windvel, N,spacing,L, sec_rot, a_wake, Nrotations, phase_diff)
TipLocation_R = rotor.TipLocation_R; 
RootLocation_R = rotor.RootLocation_R; 
TSR = rotor.TSR;        
Radius = rotor.Radius;    
NBlades = rotor.NBlades;    
[r_R] = RadialSpacing(N, TipLocation_R, RootLocation_R, spacing);

Omega = norm(windvel)*TSR/Radius;

theta_array = [0:pi/10:2*pi*Nrotations];%Omega*t, where t is the time
% Lw_D:  wake length in diameters downstream
Lw_D = max(theta_array)/Omega*norm(windvel)*(1-a_wake)/(2*Radius); % [-]

%% LLT calculations
RotorWakeSystem = vortex_system(r_R, Radius, TSR/(1-a_wake), theta_array, NBlades,sec_rot,L,phase_diff);

[InfluenceMatrix] = influence_matrix(RotorWakeSystem);
[a, aline, r_R_cp, Fnorm, Ftan, GammaNew, alpha, inflow, conv]= solveSystem(InfluenceMatrix, RotorWakeSystem, Radius, Omega, windvel,L, sec_rot,N, NBlades);

if sec_rot==1
    r_R_cp = [r_R_cp(1:N*NBlades);r_R_cp(1:N*NBlades)];
end
[CT, CP, CQ, ct, cp, cq] = CT_CPcalculations(N,Fnorm, Ftan, windvel(1), r_R, r_R_cp, Omega, Radius, NBlades, sec_rot);

%plotting_func(sec_rot,windvel,Radius, N, NBlades, Omega, a, aline, r_R_cp, ct,cp,cq, GammaNew, alpha, inflow);

results = [r_R_cp,a,aline,ct,cp,cq,GammaNew,alpha',inflow', Fnorm, Ftan];

end