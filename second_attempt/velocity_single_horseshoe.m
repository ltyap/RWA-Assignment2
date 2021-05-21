function [InducedVelocity] = velocity_single_horseshoe(localcp, localring, Gamma)
    CORE = 1e-5;
    InducedVelocity = zeros(3,1);
    Nfil = length(localring.x)-1;%number of filaments in each ring
    %start with trailing vortex along point x1 of bound vortex
    for ifil = 1:Nfil
        XV1 = [localring.x(ifil),localring.y(ifil),localring.z(ifil)];
        XV2 = [localring.x(ifil+1),localring.y(ifil+1),localring.z(ifil+1)];
        tempVelocity = v_3D_from_vortex(Gamma, XV1, XV2, localcp, CORE);        
        InducedVelocity = InducedVelocity + tempVelocity;
    end 
end

