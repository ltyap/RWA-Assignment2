function [CTrotor, CProtor, CQrotor, ct_r, cp_r, cq_r] = CT_CPcalculations(Fnorm,Ftan, Uinf, r_R, r_R_cp, Omega, Radius, NBlades)

drtemp = r_R(2:end)-r_R(1:end-1);
drtemp = repmat(drtemp,1,NBlades);
drtemp = drtemp(:);
ct_r = (drtemp.*Fnorm)/(0.5*Uinf*Uinf*pi*Radius);
cp_r = (drtemp.*Ftan.*r_R_cp.*Omega)/(0.5*(Uinf^3)*pi);
cq_r = cp_r/Omega;
CTrotor = sum(ct_r);
CProtor = sum(cp_r);
CQrotor = sum(cq_r);
end

