function [UMat, VMat, WMat] = InfluenceMatrix(RotorWakeSystem, NBlades)
    controlpoints = RotorWakeSystem.cp;
    bound = RotorWakeSystem.bound;
    trail = RotorWakeSystem.trail;
    Npan = RotorWakeSystem.NpanelsPerBlade;
    Gamma = 1;  
    Ncp = length(controlpoints.x);   
    
    UMat = zeros(Ncp,Npan);
    VMat = zeros(Ncp,Npan);
    WMat = zeros(Ncp,Npan);
    
    idx = [1:Npan+1];
    for icp = 1:Ncp
        localcp = [controlpoints.x(icp);...
                   controlpoints.y(icp);...
                   controlpoints.z(icp)];
       for iBlade = 1:NBlades
           %bound and trailing vortices on each blade
           bound_blades.x = bound.x(idx+(iBlade-1)*(Npan+1));
           bound_blades.y = bound.y(idx+(iBlade-1)*(Npan+1));
           bound_blades.z = bound.z(idx+(iBlade-1)*(Npan+1));
           trail_blades.x = trail.x(:,idx+(iBlade-1)*(Npan+1));
           trail_blades.y = trail.y(:,idx+(iBlade-1)*(Npan+1));
           trail_blades.z = trail.z(:,idx+(iBlade-1)*(Npan+1));
           for jring = 1:Npan
               bound_j.point1 = [bound_blades.x(jring),bound_blades.y(jring),bound_blades.z(jring)];
               bound_j.point2 = [bound_blades.x(jring+1),bound_blades.y(jring+1),bound_blades.z(jring+1)];
               trail_j.point1.x = trail_blades.x(:,jring);
               trail_j.point1.y = trail_blades.y(:,jring);
               trail_j.point1.z = trail_blades.z(:,jring);
               
               trail_j.point2.x = trail_blades.x(:,jring+1);
               trail_j.point2.y = trail_blades.y(:,jring+1);
               trail_j.point2.z = trail_blades.z(:,jring+1);              
               
               InducedVelocity = velocity_single_horseshoe(localcp, bound_j, trail_j, Gamma);
               UMat(icp,jring+(iBlade-1)*Npan) = InducedVelocity(1);
               VMat(icp,jring+(iBlade-1)*Npan) = InducedVelocity(2);
               WMat(icp,jring+(iBlade-1)*Npan) = InducedVelocity(3);
           end
       end
    end

end
