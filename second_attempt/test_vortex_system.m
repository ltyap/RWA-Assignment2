function [results]=vortex_system(r_R, Radius, tipspeedratio, theta_array, NBlades,sec_rot)
[chord_distribution, twist_distribution] = BladeGeometry(r_R);
chord_cp = (chord_distribution(1:end-1)+chord_distribution(2:end))/2;
N = length(r_R)-1;  % number of panels/control points/bound vortices
trail.x = [];
trail.y = [];
trail.z = [];
for n=1:NBlades % loop over blades
    blade_azim = 2*pi/NBlades*(n-1);
    
    % determine control point locations for no-normal flow condition - at middle of segments
    r_cp = (r_R(1:end-1)+r_R(2:end))/2*Radius;%[m]
    [~,twist] = BladeGeometry(r_cp/Radius); % twist at control points [deg]
    twist = deg2rad(twist);
    cp.coordinates(1,(n-1)*N+1:n*N) = zeros(N,1);%0.5*chord_cp';    % x-coordinate
    cp.coordinates(2,(n-1)*N+1:n*N) = r_cp.*cos(blade_azim);    % y-coordinate
    cp.coordinates(3,(n-1)*N+1:n*N) = r_cp.*sin(blade_azim);    % z-coordinate
    cp.radius((n-1)*N+1:n*N) = r_cp;
    cp.chord((n-1)*N+1:n*N) = chord_cp;
    
    cp.normal(:,(n-1)*N+1:n*N) = [cos(twist);...
        0*cos(blade_azim) - -sin(twist)*sin(blade_azim);...
        0*sin(blade_azim) + -sin(twist)*cos(blade_azim)];
    cp.tangential(:,(n-1)*N+1:n*N) = [-sin(twist); ...
        0*cos(blade_azim) - -cos(twist)*sin(blade_azim);...
        0*sin(blade_azim) + -cos(twist)*cos(blade_azim)];
    
    % bound vortices
    bound.x((n-1)*(N+1)+1:n*(N+1)) = zeros(N+1,1); % x-coordinate
    bound.y((n-1)*(N+1)+1:n*(N+1)) = r_R*Radius*cos(blade_azim);    % y-coordinate
    bound.z((n-1)*(N+1)+1:n*(N+1)) = r_R*Radius*sin(blade_azim);    % z-coordinate
    
    bound.centrepoint(1,(n-1)*N+1:n*N) = zeros(N,1);
    bound.centrepoint(2,(n-1)*N+1:n*N) = r_cp*cos(blade_azim);
    bound.centrepoint(3,(n-1)*N+1:n*N) = r_cp*sin(blade_azim);
    
    bound.radialcentrepoint((n-1)*N+1:n*N) = sqrt(dot(bound.centrepoint(:,(n-1)*N+1:n*N),...
        bound.centrepoint(:,(n-1)*N+1:n*N)));
    
    % trailing vortices
    for i=1:N+1 % loop over bound vortices
        r = r_R(i);
        c = chord_distribution(i);
        beta = -deg2rad(twist_distribution(i));
        
        % point 1 of the 1st pair of trailing vortices defined by bound
        % vortices, so only need to define 2nd point
        temp1.x(i) = bound.x((n-1)*(N+1)+i)+c*sin(beta)
        temp1.y(i) = bound.y((n-1)*(N+1)+i)
        temp1.z(i) = bound.z((n-1)*(N+1)+i)+c*cos(beta)
        
        %Rotate them, depending on blade
        temp1.y(i) = temp1.y(i)*cos(blade_azim) - temp1.z(i)*sin(blade_azim);
        temp1.z(i) = temp1.y(i)*sin(blade_azim) + temp1.z(i)*cos(blade_azim);       
        
        
        trail.x = [trail.x,temp1.x(i)];
        trail.y = [trail.y,temp1.y(i)];
        trail.z = [trail.z,temp1.z(i)];
        
        for j = 1:length(theta_array)-1
            dx = (theta_array(j+1)-theta_array(j))/tipspeedratio*Radius;
            dy = r_R(i)*Radius*(cos(-theta_array(j+1)+blade_azim)-cos(-theta_array(j)+blade_azim));%
            dz = r_R(i)*Radius*(sin(-theta_array(j+1)+blade_azim)-sin(-theta_array(j)+blade_azim));%
            temp1.x = trail.x(end)+dx;
            temp1.y = trail.y(end)+dy;
            temp1.z = trail.z(end)+dz;
            trail.x = [trail.x, temp1.x];
            trail.y = [trail.y, temp1.y];
            trail.z = [trail.z, temp1.z];
        end
        
        
    end
end

trail.x = reshape(trail.x,length(theta_array),NBlades*length(r_R));
trail.y = reshape(trail.y,length(theta_array),NBlades*length(r_R));
trail.z = reshape(trail.z,length(theta_array),NBlades*length(r_R));


% plot system - for checking
    hold on
% scatter3(cp.coordinates(1,:), cp.coordinates(2,:), cp.coordinates(3,:), 10, 'filled', 'r')
scatter3(bound.centrepoint(1,:), bound.centrepoint(2,:), bound.centrepoint(3,:),'r*');
xlabel("x")
ylabel("y")
zlabel("z")
xlim([-20 20])
for n=1:NBlades
    plot3(bound.x((n-1)*(N+1)+1:n*(N+1)), bound.y((n-1)*(N+1)+1:n*(N+1)), bound.z((n-1)*(N+1)+1:n*(N+1)), '-k+')
    for i=1:N+1
        plot3([trail.x(1,(n-1)*(N+1)+i), bound.x((n-1)*(N+1)+i)], ...
            [trail.y(1,(n-1)*(N+1)+i), bound.y((n-1)*(N+1)+i)], [trail.z(1,(n-1)*(N+1)+i), bound.z((n-1)*(N+1)+i)], '--bo')
        plot3(trail.x(:,(n-1)*(N+1)+i), trail.y(:,(n-1)*(N+1)+i), trail.z(:,(n-1)*(N+1)+i), '--r')
    end
end
%     hold off


trail.x = [bound.x;trail.x];
trail.y = [bound.y;trail.y];
trail.z = [bound.z;trail.z];

trail.x = flip(trail.x);
trail.y = flip(trail.y);
trail.z = flip(trail.z);
ring.x = [];
ring.y = [];
ring.z = [];
idxfilaments = [1:N+1];
for n = 1:NBlades
    filaments.x = trail.x(:,idxfilaments+(n-1)*(N+1));
    filaments.y = trail.y(:,idxfilaments+(n-1)*(N+1));
    filaments.z = trail.z(:,idxfilaments+(n-1)*(N+1));
    for i = 1:N
        temp2.x(:,i) = [filaments.x(:,i);flip(filaments.x(:,i+1))];
        temp2.y(:,i) = [filaments.y(:,i);flip(filaments.y(:,i+1))];
        temp2.z(:,i) = [filaments.z(:,i);flip(filaments.z(:,i+1))];
    end
    ring.x = [ring.x, temp2.x];
    ring.y = [ring.y, temp2.y];
    ring.z = [ring.z, temp2.z];
end

if sec_rot==1 % 2nd rotor
    SeparationDist = 2*(2*Radius);
    phase = deg2rad(80);
    cosphase = cos(phase);
    sinephase = sin(phase);
    
    transformed(1,:) = cp.coordinates(1,:);
    transformed(2,:) = (cp.coordinates(2,:))*cosphase - cp.coordinates(3,:)*sinephase +SeparationDist;
    transformed(3,:) = (cp.coordinates(2,:))*sinephase + cp.coordinates(3,:)*cosphase;
    scatter3(transformed(1,:), transformed(2,:), transformed(3,:), 10, 'filled', 'r')
    cp.coordinates = [cp.coordinates, transformed];
    
    clear transformed
    transformed(1,:) = bound.centrepoint(1,:);
    transformed(2,:) = (bound.centrepoint(2,:))*cosphase - bound.centrepoint(3,:)*sinephase +SeparationDist;
    transformed(3,:) = (bound.centrepoint(2,:))*sinephase + bound.centrepoint(3,:)*cosphase;
    scatter3(transformed(1,:),transformed(2,:), transformed(3,:),'r*');
    bound.centrepoint = [bound.centrepoint, transformed];
    
    
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

cp.Totalcp = N*NBlades;
results.ring = ring;
results.NBlades = NBlades;
results.cp = cp;
results.bound = bound;
results.trail = trail;
results.NpanelsPerBlade = N;
end