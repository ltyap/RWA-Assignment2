function [InfluenceMatrix] = InfluenceMatrix(RotorWakeSystem, NBlades)
%     controlpoints = RotorWakeSystem.cp;
    bound = RotorWakeSystem.bound;
%     trail = RotorWakeSystem.trail;
    ring = RotorWakeSystem.ring;
    Npan = RotorWakeSystem.NpanelsPerBlade;
    Gamma = 1;  
    Ncp = RotorWakeSystem.cp.Totalcp;   
    
    UMat = zeros(Ncp,Npan*NBlades);
    VMat = zeros(Ncp,Npan*NBlades);
    WMat = zeros(Ncp,Npan*NBlades);
    
    idx = [1:Npan];
    for icp = 1:Ncp
        localcp = bound.centrepoint(:,icp);
        for iBlade = 1:NBlades
            RingOnBlade.x = ring.x(:,idx+(iBlade-1)*Npan);
            RingOnBlade.y = ring.y(:,idx+(iBlade-1)*Npan);
            RingOnBlade.z = ring.z(:,idx+(iBlade-1)*Npan);
            for jring = 1:Npan
                localring.x = RingOnBlade.x(:,jring);
                localring.y = RingOnBlade.y(:,jring);
                localring.z = RingOnBlade.z(:,jring);                          
               InducedVelocity = velocity_single_horseshoe(localcp, localring, Gamma);
               UMat(icp,jring+(iBlade-1)*Npan) = InducedVelocity(1);
               VMat(icp,jring+(iBlade-1)*Npan) = InducedVelocity(2);
               WMat(icp,jring+(iBlade-1)*Npan) = InducedVelocity(3);
            end
        end
    end
    
    InfluenceMatrix.U = UMat;
    InfluenceMatrix.V = VMat;
    InfluenceMatrix.W = WMat;
%% Plotting influence matrices - for checking
%     figure;
%     imagesc(UMat);
%     colorbar();
%     title('U matrix');
%     figure;
%     imagesc(VMat);
%     title('Vmatrix');
%     figure;
%     imagesc(WMat);
%     title('W Matrix');     
end