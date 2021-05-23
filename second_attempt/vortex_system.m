function [results]=vortex_system(r_R, Radius, tipspeedratio, theta_array, NBlades,sec_rot)
[c, twist_distribution] = BladeGeometry(r_R);
twist = -deg2rad(twist_distribution);
cp_chord = (c(1:end-1)+c(2:end))/2;
r_R_cp = 0.5*(r_R(1:end-1) + r_R(2:end));

[~,cp_twist] = BladeGeometry(r_R_cp);
cp_twist = deg2rad(cp_twist);

theta_array = theta_array(:);
N = length(r_R)-1;  % number of panels/control points/bound vortices

%start with one blade, then rotate to get points on the rest of the blades
%% control points at centre of bound vortex/0.75 chord
bound.cpcoord(1,:) = zeros(1,N);%x
bound.cpcoord(2,:) = 0.5*(r_R(2:end)+r_R(1:end-1))*Radius;%y
bound.cpcoord(3,:) = zeros(1,N);%z

%checking
scatter3(bound.cpcoord(1,:),bound.cpcoord(2,:), bound.cpcoord(3,:), 'bo');
xlabel('x');
ylabel('y');
zlabel('z');
hold on
%% normal and tangential vector of control point (for no normal flow condition)
bound.normal = [cos(cp_twist);zeros(1,N);-sin(cp_twist)];
bound.tangential = [-sin(cp_twist);zeros(1,N);-cos(cp_twist)];
%% bound vortices
bound.x = zeros(1,N+1);
bound.y = r_R*Radius;
bound.z = zeros(1,N+1);

for i = 1:N
    plot3([bound.x(i),bound.x(i+1)],...
        [bound.y(i),bound.y(i+1)],...
        [bound.z(i),bound.z(i+1)],'k-x');
end
%% trailing vortices
%first pair of trailing vortices
trail.x = zeros(length(theta_array)+1,N+1);
trail.y = zeros(length(theta_array)+1,N+1);
trail.z = zeros(length(theta_array)+1,N+1);

trail.x(1,:) = bound.x;
trail.y(1,:) = bound.y;
trail.z(1,:) = bound.z;

trail.x(2,:) = trail.x(1,:) - c.*sin(twist);
trail.y(2,:) = trail.y(1,:);
trail.z(2,:) = trail.z(1,:) - c.*cos(twist);

for i = 1:N+1
    plot3([trail.x(1,i), trail.x(2,i)],[trail.y(1,i), trail.y(2,i)],...
        [trail.z(1,i), trail.z(2,i)],'k--+');
end

%trailing vortices in wake
for i = 1:N+1
    for j = 1:length(theta_array)-1
        dx = (theta_array(j+1)-theta_array(j))/tipspeedratio*Radius;
        dy = r_R(i)*Radius*(cos(-theta_array(j+1))-cos(-theta_array(j)));%
        dz = r_R(i)*Radius*(sin(-theta_array(j+1))-sin(-theta_array(j)));%
        trail.x(j+2,i) = trail.x(j+1,i) + dx;
        trail.y(j+2,i) = trail.y(j+1,i) + dy; 
        trail.z(j+2,i) = trail.z(j+1,i) + dz; 
    end
end

Nfil = length(trail.x(:,1));

for i = 1:N+1
    plot3(trail.x(:,i), trail.y(1:end,i), trail.z(1:end,i),'r--');
end
%% Ring
for i = 1:N
    ring.x(:,i) = [flip(trail.x(:,i));trail.x(:,i+1)];
    ring.y(:,i) = [flip(trail.y(:,i));trail.y(:,i+1)];
    ring.z(:,i) = [flip(trail.z(:,i));trail.z(:,i+1)];
end

%% Rotation
bladeAzim = linspace(0,2*pi,NBlades+1);
bladeAzim = bladeAzim(2:end-1);

for krot = 1:length(bladeAzim)
    for i = 1:N
        newboundcpcoord(:,i) = Rotate(bound.cpcoord(:,i),bladeAzim(krot));
        newnormal(:,i) = Rotate(bound.normal(:,i),bladeAzim(krot));
        newtangential(:,i) = Rotate(bound.tangential(:,i),bladeAzim(krot));
    end
    bound.normal = [bound.normal, newnormal];
    bound.tangential = [bound.tangential, newtangential];
    bound.cpcoord = [bound.cpcoord, newboundcpcoord];
    scatter3(newboundcpcoord(1,:),newboundcpcoord(2,:),newboundcpcoord(3,:),'bo');
end
for krot = 1:length(bladeAzim)
    rotatedrings_x = ring.x;
    rotatedrings_y = ring.y*cos(bladeAzim(krot)) - ring.z*sin(bladeAzim(krot));% y positions
    rotatedrings_z = ring.y*sin(bladeAzim(krot)) + ring.z*cos(bladeAzim(krot));% z positions
    for i = 1:N   
        plot3(rotatedrings_x(1:Nfil,i),...
            rotatedrings_y(1:Nfil,i),...
            rotatedrings_z(1:Nfil,i),'r--');
        plot3(rotatedrings_x(Nfil:Nfil+1,i),...
            rotatedrings_y(Nfil:Nfil+1,i),...
            rotatedrings_z(Nfil:Nfil+1,i),'k-x');
        plot3(rotatedrings_x(Nfil+1:end,i),...
            rotatedrings_y(Nfil+1:end,i),...
            rotatedrings_z(Nfil+1:end,i),'r--');
    end
    ring.x = [ring.x, rotatedrings_x];
    ring.y = [ring.y, rotatedrings_y];
    ring.z = [ring.z, rotatedrings_z];
end
if sec_rot == 1
    bound.Totalcp = N*NBlades*2;
elseif sec_rot == 0
    bound.Totalcp = N*NBlades;
else
    disp('Invalid value of sec_rot, proceeding to use single rotor only');
    bound.Totalcp = N*NBlades;
end
 
if sec_rot==1 % 2nd rotor
    SeparationDist = 2*(2*Radius);
    phase = deg2rad(80);
    cosphase = cos(phase);
    sinephase = sin(phase);
    
    
    transformed(1,:) = bound.normal(1,:);
    transformed(2,:) = (bound.normal(2,:))*cosphase - bound.normal(3,:)*sinephase ;
    transformed(3,:) = (bound.normal(2,:))*sinephase + bound.normal(3,:)*cosphase;
    bound.normal = [bound.normal, transformed];
    clear transformed
    
    transformed(1,:) = bound.tangential(1,:);
    transformed(2,:) = (bound.tangential(2,:))*cosphase - bound.tangential(3,:)*sinephase;
    transformed(3,:) = (bound.tangential(2,:))*sinephase + bound.tangential(3,:)*cosphase;
    bound.tangential = [bound.tangential, transformed];
    clear transformed
   
    clear transformed
    transformed(1,:) = bound.cpcoord(1,:);
    transformed(2,:) = (bound.cpcoord(2,:))*cosphase - bound.cpcoord(3,:)*sinephase +SeparationDist;
    transformed(3,:) = (bound.cpcoord(2,:))*sinephase + bound.cpcoord(3,:)*cosphase;
    scatter3(transformed(1,:),transformed(2,:), transformed(3,:),'o');
    bound.cpcoord = [bound.cpcoord, transformed];
    
    
    clear transformed
    transformed(:,:,1) = ring.x;
    transformed(:,:,2) = ring.y*cosphase - ring.z*sinephase + SeparationDist;
    transformed(:,:,3) = ring.y*sinephase + ring.z*cosphase;
    for n=1:NBlades
        for i=1:N
            plot3(transformed(:,(n-1)*N+i,1),...
                transformed(:,(n-1)*N+i,2),...
                transformed(:,(n-1)*N+i,3), 'b');
        end
    end
    ring.x = [ring.x, transformed(:,:,1)];
    ring.y = [ring.y, transformed(:,:,2)];
    ring.z = [ring.z, transformed(:,:,3)];
end
 
bound.cp_radialpos = sqrt(dot(bound.cpcoord,bound.cpcoord,1));

% checkvectors = dot(bound.normal, bound.tangential,1);


results.ring = ring;
results.NBlades = NBlades;%blade on each rotor for now
results.bound = bound;
results.NpanelsPerBlade = N;
end