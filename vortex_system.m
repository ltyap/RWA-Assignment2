function [cp, bound, trail]=vortex_system(Omega, NBlades, theta_array,N, r_R, Radius, twist_distribution, chord_distribution, Lw_D, Uwake)
ts = 0.2; % maybe come up with a function for determining this
t = [0:ts:Lw_D*2*Radius/Uwake];

for n=1:NBlades % loop over blades
    blade_azim = theta_array(n);
    
    % determine control point locations - at middle of segments
    r = (r_R(1:end-1)+r_R(2:end))/2*Radius;
    cp.x((n-1)*N+1:n*N) = zeros(N,1); % x-coordinate
    cp.y((n-1)*N+1:n*N) = r*cos(2*pi-blade_azim);     % y-coordinate
    cp.z((n-1)*N+1:n*N) = r*sin(2*pi-blade_azim);     % z-coordinate
    cp.twist((n-1)*N+1:n*N) = interp1(r_R, twist_distribution, r/Radius); % twist at control points
    
    % bound vortices
    bound.x((n-1)*(N+1)+1:n*(N+1)) = zeros(N+1,1); % x-coordinate
    bound.y((n-1)*(N+1)+1:n*(N+1)) = r_R*Radius*cos(2*pi-blade_azim);     % y-coordinate
    bound.z((n-1)*(N+1)+1:n*(N+1)) = r_R*Radius*sin(2*pi-blade_azim);     % z-coordinate
    
    %determine trailing vortices
    for i=1:N+1 % loop over bound vortices
        r = r_R(i);
        c = chord_distribution(i);
        beta = -deg2rad(twist_distribution(i));
        
        % 1st trailing vortex, not sure about this
        trailx((n-1)*(N+1)+i) = bound.x((n-1)*(N+1)+i)+c*cos(beta);
        traily((n-1)*(N+1)+i) = bound.y((n-1)*(N+1)+i)-c*sin(beta)*sin(blade_azim);
        trailz((n-1)*(N+1)+i) = bound.z((n-1)*(N+1)+i)+c*sin(beta)*cos(blade_azim);
        
        x = t*Uwake;
        y = r_R(i)*Radius*sin(Omega*t+pi/2+blade_azim);
        z = r_R(i)*Radius*cos(Omega*t+pi/2+blade_azim);
        trail.x((n-1)*(N+1)+i,:) = x+trailx((n-1)*(N+1)+i);
        trail.y((n-1)*(N+1)+i,:) = y-c*sin(beta)*sin(blade_azim); %[trail_y1, trail_y1+y];
        trail.z((n-1)*(N+1)+i,:) = z+c*sin(beta)*cos(blade_azim); %[trail_z1, trail_z1+z];
    end
end

% plot system - for checking
hold on
scatter3(cp.x, cp.y, cp.z, 10, 'filled', 'r')
for n=1:NBlades
    plot3(bound.x((n-1)*(N+1)+1:n*(N+1)), bound.y((n-1)*(N+1)+1:n*(N+1)), bound.z((n-1)*(N+1)+1:n*(N+1)), '-k+')
    for i=1:N+1
        plot3([trail.x((n-1)*(N+1)+i), bound.x((n-1)*(N+1)+i)], [trail.y((n-1)*(N+1)+i), bound.y((n-1)*(N+1)+i)], [trail.z((n-1)*(N+1)+i), bound.z((n-1)*(N+1)+i)], '-k+')
        plot3(trail.x((n-1)*(N+1)+i,:), trail.y((n-1)*(N+1)+i,:), trail.z((n-1)*(N+1)+i,:), '-b+')
    end
end
% for i=1:length(trail2.x(:,1))
%     plot3(trail2.x(i,:), trail2.y(i,:), trail2.z(i,:), '-k+')
% end
% axis equal
hold off
xlabel("x")
ylabel("y")
zlabel("z")
end