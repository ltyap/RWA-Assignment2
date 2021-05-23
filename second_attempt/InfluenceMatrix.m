function [InfluenceMatrix] = InfluenceMatrix(RotorWakeSystem)
bound = RotorWakeSystem.bound;
ring = RotorWakeSystem.ring;
Npan = RotorWakeSystem.NpanelsPerBlade;
Gamma = 1;
Ncp = bound.Totalcp;
NBlades = RotorWakeSystem.NBlades;

UMat = zeros(Ncp,Ncp);
VMat = zeros(Ncp,Ncp);
WMat = zeros(Ncp,Ncp);

idx = [1:Npan];
    for icp = 1:Ncp
        localcp = bound.cpcoord(:,icp);
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

% for icp = 1:Ncp
%     localcp = bound.cpcoord(:,icp);
%     for iRotor = 1:sec_rot+1
%         for iBlade = 1:NBlades
%             RingOnBlade.x = ring.x(:,idx+(iBlade-1)*Npan+(iRotor-1)*Npan*NBlades);
%             RingOnBlade.y = ring.y(:,idx+(iBlade-1)*Npan+(iRotor-1)*Npan*NBlades);
%             RingOnBlade.z = ring.z(:,idx+(iBlade-1)*Npan+(iRotor-1)*Npan*NBlades);
%             for jring = 1:Npan
%                 localring.x = RingOnBlade.x(:,jring);
%                 localring.y = RingOnBlade.y(:,jring);
%                 localring.z = RingOnBlade.z(:,jring);
%                 InducedVelocity = velocity_single_horseshoe(localcp, localring, Gamma);
%                 UMat(icp,jring+(iBlade-1)*Npan+(iRotor-1)*Npan*NBlades) = InducedVelocity(1);
%                 VMat(icp,jring+(iBlade-1)*Npan+(iRotor-1)*Npan*NBlades) = InducedVelocity(2);
%                 WMat(icp,jring+(iBlade-1)*Npan+(iRotor-1)*Npan*NBlades) = InducedVelocity(3);
%             end
%         end
%     end
% end

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