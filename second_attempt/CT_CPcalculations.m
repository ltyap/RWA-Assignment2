function [CTrotor, CProtor, CQrotor] = CT_CPcalculations(Fnorm,Ftan, Uinf, r_R, r_R_cp, Omega, Radius, NBlades)

drtemp = r_R(2:end)-r_R(1:end-1);
drtemp = repmat(drtemp,1,NBlades);
drtemp = drtemp(:);
CTtemp = (drtemp.*Fnorm)/(0.5*Uinf*Uinf*pi*Radius);
CPtemp = (drtemp.*Ftan.*r_R_cp.*Omega)/(0.5*(Uinf^3)*pi);
CQtemp = CPtemp/Omega;
CTrotor = sum(CTtemp);
CProtor = sum(CPtemp);
CQrotor = sum(CQtemp);
end

