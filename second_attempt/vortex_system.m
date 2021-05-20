function [results]=vortex_system(r_R, Radius, tipspeedratio, theta_array, NBlades)
    [chord_distribution, twist_distribution] = BladeGeometry(r_R);
    N = length(r_R)-1;
    trail.x = [];
    trail.y = [];
    trail.z = [];
%     dtheta = theta_array(2:end)-theta_array(1:end-1);
    for n=1:NBlades % loop over blades
        blade_azim = 2*pi/NBlades*(n-1);

        % determine control point locations - at middle of segments
        r_cp = (r_R(1:end-1)+r_R(2:end))/2*Radius;%[m]
        chord_cp = (chord_distribution(1:end-1)+chord_distribution(2:end))/2;
        [~,temp] = BladeGeometry(r_cp/Radius); % twist at control points [deg]
        cp.twist((n-1)*N+1:n*N) = temp;
        cp.x((n-1)*N+1:n*N) = zeros(N,1)+(0.75*chord_cp'); % x-coordinate
        cp.y((n-1)*N+1:n*N) = r_cp*cos(blade_azim);     % y-coordinate
        cp.z((n-1)*N+1:n*N) = r_cp.*cosd(temp)*sin(blade_azim);     % z-coordinate

        % bound vortices
        bound.x((n-1)*(N+1)+1:n*(N+1)) = zeros(N+1,1)+0.25*chord_distribution'; % x-coordinate
        bound.y((n-1)*(N+1)+1:n*(N+1)) = r_R*Radius*cos(blade_azim);     % y-coordinate
        bound.z((n-1)*(N+1)+1:n*(N+1)) = r_R*Radius*sin(blade_azim);     % z-coordinate
        %determine trailing vortices
        for i=1:N+1 % loop over bound vortices
            r = r_R(i);
            c = chord_distribution(i);
            beta = -deg2rad(twist_distribution(i));

            % 1st trailing vortex, not sure about this
            temp1.x(i) = bound.x((n-1)*(N+1)+i)+0.75*c*cos(beta);%(n-1)*(N+1)+
            temp1.y(i) = bound.y((n-1)*(N+1)+i)-0.75*c*sin(beta)*sin(blade_azim);
            temp1.z(i) = bound.z((n-1)*(N+1)+i)+0.75*c*sin(beta)*cos(blade_azim);
                        
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
    scatter3(cp.x, cp.y, cp.z, 10, 'filled', 'r')
    for n=1:NBlades
        plot3(bound.x((n-1)*(N+1)+1:n*(N+1)), bound.y((n-1)*(N+1)+1:n*(N+1)), bound.z((n-1)*(N+1)+1:n*(N+1)), '-k+')
        for i=1:N+1
            plot3([trail.x(1,(n-1)*(N+1)+i), bound.x((n-1)*(N+1)+i)], ...
                [trail.y(1,(n-1)*(N+1)+i), bound.y((n-1)*(N+1)+i)], [trail.z(1,(n-1)*(N+1)+i), bound.z((n-1)*(N+1)+i)], '--bo')
            plot3(trail.x(:,(n-1)*(N+1)+i), trail.y(:,(n-1)*(N+1)+i), trail.z(:,(n-1)*(N+1)+i), '-r')
        end
    end
    hold off
    xlabel("x")
    ylabel("y")
    zlabel("z")
    
    results.cp = cp;
    results.bound = bound;
    results.trail = trail;
end