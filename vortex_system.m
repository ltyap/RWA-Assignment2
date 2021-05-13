function [cp, bound, trail]=vortex_system(NBlades, theta_array,N, r_R, Radius, twist_distribution, chord_distribution)

for n=1:NBlades
    blade_azim = theta_array(n);
    
    % determine control point locations - at middle of segments
    r = (r_R(1:end-1)+r_R(2:end))/2*Radius;
    cp.x((n-1)*N+1:n*N) = zeros(N,1); % x-coordinate
    cp.y((n-1)*N+1:n*N) = r*cos(blade_azim);     % y-coordinate
    cp.z((n-1)*N+1:n*N) = r*sin(blade_azim)     % z-coordinate
    cp.twist((n-1)*N+1:n*N) = interp1(r_R, twist_distribution, r/Radius); % twist at control points
    
    % bound vortices
    bound.x((n-1)*(N+1)+1:n*(N+1)) = zeros(N+1,1); % x-coordinate
    bound.y((n-1)*(N+1)+1:n*(N+1)) = r_R*Radius*cos(blade_azim);     % y-coordinate
    bound.z((n-1)*(N+1)+1:n*(N+1)) = r_R*Radius*sin(blade_azim)     % z-coordinate
    
    %determine trailing vortices
    for i=1:N+1 % loop over bound vortices
        r = r_R(i);
        c = chord_distribution(i);
        beta = -deg2rad(twist_distribution(i));
        
        % baseline
        trail.x((n-1)*(N+1)+i) = bound.x((n-1)*(N+1)+i)+c*cos(beta);
        trail.y((n-1)*(N+1)+i) = bound.y((n-1)*(N+1)+i)-c*sin(beta)*sin(blade_azim);
        trail.z((n-1)*(N+1)+i) = bound.z((n-1)*(N+1)+i)+c*sin(beta)*cos(blade_azim);
        
        % trailing vortices
        %                 trail_x1 = 0 + chord_distribution(i)*sin(deg2rad(-twist_distribution(i)));
        %                 trail_y1 = chord_distribution(i)/4 + chord_distribution(i)*cos(deg2rad(-twist_distribution(i)));
        %                 trail_z1 = r_R(i)*Radius + 0;
        %
        
        
        %ts_func = @(ts) sqrt(ts^2*Uwake^2+r_R(i)^2*Radius^2*sin(Omega*ts)^2+r_R(i)^2*Radius^2*cos(Omega*ts)^2)-chord_distribution(i) ;
        %ts = fsolve(ts_func,0.1);  % step so that we have wake discretized into chord lengths
        %                 ts = 0.2; % maybe come up with a function for determining this
        %                 t = [0:ts:Lw_D*2*Radius/Uwake];
        %                 x = t*Uwake;
        %                 y = r_R(i)*Radius*sin(Omega*t);
        %                 z = r_R(i)*Radius*cos(Omega*t);
        %                 trail.x(i,:) = [trail_x1, trail_x1+x];
        %                 trail.y(i,:) = [trail_y1, trail_y1+y]; %[trail_y1, trail_y1+y];
        %                 trail.z(i,:) = [trail_z1, z]; %[trail_z1, trail_z1+z];
    end
    
    
end

% plot system - for checking
hold on
scatter3(cp.x, cp.y, cp.z, 10, 'filled', 'r')
for n=1:NBlades
    plot3(bound.x((n-1)*(N+1)+1:n*(N+1)), bound.y((n-1)*(N+1)+1:n*(N+1)), bound.z((n-1)*(N+1)+1:n*(N+1)), '-k+')
    for i=1:N+1
        plot3([trail.x((n-1)*(N+1)+i), bound.x((n-1)*(N+1)+i)], [trail.y((n-1)*(N+1)+i), bound.y((n-1)*(N+1)+i)], [trail.z((n-1)*(N+1)+i), bound.z((n-1)*(N+1)+i)], '-k+')
    end
end

hold off
xlabel("x")
ylabel("y")
zlabel("z")
end