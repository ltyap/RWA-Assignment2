function [CTrotor, CProtor, CQrotor, ct_r, cp_r, cq_r] = CT_CPcalculations(N,Fnorm,Ftan, Uinf, r_R, r_R_cp, Omega, Radius, NBlades, sec_rot)

drtemp = r_R(2:end)-r_R(1:end-1);
drtemp = repmat(drtemp,1,NBlades);
drtemp = drtemp(:);
if sec_rot==1
    drtemp = [drtemp; drtemp];
end
ct_r = (drtemp.*Fnorm*NBlades)/(0.5*Uinf^2*pi*Radius^2);
cp_r = (drtemp.*Ftan.*r_R_cp.*Omega*Radius*NBlades)/(0.5*(Uinf^3)*pi*Radius^2);
cq_r = cp_r/Omega;

if sec_rot==1
    CTrotor = [sum(ct_r(1:N*NBlades)), sum(ct_r(N*NBlades+1:end))];
    CProtor = [sum(cp_r(1:N*NBlades)), sum(cp_r(N*NBlades+1:end))];
    CQrotor = [sum(cq_r(1:N*NBlades)), sum(cq_r(N*NBlades+1:end))];
else
    CTrotor = sum(ct_r);
    CProtor = sum(cp_r);
    CQrotor = sum(cq_r);
end
end

