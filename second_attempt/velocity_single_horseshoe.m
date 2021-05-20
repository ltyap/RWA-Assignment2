function [InducedVelocity] = velocity_single_horseshoe(localcp, bound_j, trail_j, Gamma)
    
    InducedVelocity = zeros(3,1);
    Nfil = length(trail_j.point1.x)-1;
    %start with trailing vortex at point x1 of bound vortex
    for ifil = 1:Nfil
        XV1 = [trail_j.point1.x(ifil),trail_j.point1.y(ifil),trail_j.point1.z(ifil)];
        XV2 = [trail_j.point1.x(ifil+1),trail_j.point1.y(ifil+1),trail_j.point1.z(ifil+1)];
        tempVelocity = v_3D_from_vortex(Gamma, XV1, XV2, localcp);        
        InducedVelocity(1) = InducedVelocity(1) + tempVelocity(1);
        InducedVelocity(2) = InducedVelocity(2) + tempVelocity(2);
        InducedVelocity(3) = InducedVelocity(3) + tempVelocity(3);
    end
    %Trailing vortex at point x2 of bound vortex
    for ifil = 1:Nfil
        XV1 = [trail_j.point2.x(ifil),trail_j.point2.y(ifil),trail_j.point2.z(ifil)];
        XV2 = [trail_j.point2.x(ifil+1),trail_j.point2.y(ifil+1),trail_j.point2.z(ifil+1)];
        tempVelocity = v_3D_from_vortex(Gamma, XV1, XV2, localcp);        
        InducedVelocity(1) = InducedVelocity(1) + tempVelocity(1);
        InducedVelocity(2) = InducedVelocity(2) + tempVelocity(2);
        InducedVelocity(3) = InducedVelocity(3) + tempVelocity(3);
    end
    %Bound vortex
    XV1 = bound_j.point1;
    XV2 = bound_j.point2;
    tempVelocity = v_3D_from_vortex(Gamma, XV1, XV2, localcp);        
    InducedVelocity(1) = InducedVelocity(1) + tempVelocity(1);
    InducedVelocity(2) = InducedVelocity(2) + tempVelocity(2);
    InducedVelocity(3) = InducedVelocity(3) + tempVelocity(3);
    

end

