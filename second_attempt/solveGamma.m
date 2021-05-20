% new iterative loop
function [a, aline, r_R, Fnorm, Ftan, Gamma_temp] = solveGamma(Uinf, N, Radius, cp, Influence_u, Influence_v, Influence_w, Omega, polar_alpha, polar_cl, polar_cd, chord_distribution)
GammaNew = ones(length(cp.x),1); % initial guess, defined in control points
U_inf = [Uinf,0,0];

Niterations = 1000;
errorlimit = 0.0001;
error = 1.0;

% Niterations = 10000;    % maximum number of performed iterations
% Erroriterations = 0.0000001; % error tolerance for the iterative process

for i=1:Niterations
    Gamma = GammaNew;%???? % update bound circulation
    
    % calculate velocity, circulation and loads at the controlpoints
    for icp= 1:length(cp.x)% loop over all control points
        local_cp = [cp.x(icp), cp.y(icp), cp.z(icp)];
        % determine radial position of the controlpoint;
        r = sqrt(dot(local_cp, local_cp));
        % QUICK FIX - rewrite!
        if icp<=N
            local_c = (chord_distribution(icp)+chord_distribution(icp+1))/2; % local chord
        elseif icp<=2*N
            local_c = (chord_distribution(icp-N)+chord_distribution(icp+1-N))/2; % local chord
        else
            local_c = (chord_distribution(icp-2*N)+chord_distribution(icp+1-2*N))/2; % local chord
        end
        
        u=0; % initialize velocity at control point
        v=0;
        w=0;
        
        % multiply icp line of Matrix with vector of circulation Gamma to calculate velocity at controlpoint
        for j = 1:length(Influence_u(1,:)) % loop over all vortex rings
            u = u + Influence_u(icp,j)*Gamma(j); % axial component of velocity
            v = v + Influence_v(icp,j)*Gamma(j); % y-component of velocity
            w = w + Influence_w(icp,j)*Gamma(j); % z-component of velocity
        end
        
        %         Uaxial = Uinf+u;
        %         Utan = Omega*r+v;%dot([Uinf+u, v, w],[1,0,0]);
         Urot = cross([-Omega,0,0],local_cp);
         Urel = [Uinf+u+Urot(1), v+Urot(2), w+Urot(3)];
         azimdir = cross([-1/r, 0, 0]  , local_cp); % rotational direction
         Utan = dot(azimdir, Urel); % azimuthal direction
         Uaxial =  dot([1, 0, 0] , Urel); % axial velocity
%         Uaxial = Uinf+u;
%         Utan = Omega*r + dot([Uinf+u, v, w],azimdir);
        
        [fnorm , ftan, gamma, ~, ~] = loadBladeElement(r/Radius, local_c, cp.twist(icp), polar_alpha, polar_cl, polar_cd, Uaxial, Utan);
        
        % new point of new estimate of circulation for the blade section
        GammaNew(icp) = gamma;
        
        % update output vector
        a(icp) = 1-Uaxial/Uinf;
        aline(icp) = Utan/(r*Omega)-1;
        r_R(icp) = r/Radius ;
        Fnorm(icp) = fnorm ;
        Ftan(icp) = ftan ;
        Gamma_temp(icp) = gamma;
    end % end loop control points
    
    % check convergence of solution
    ref_error = max(abs(GammaNew));
    ref_error = max(ref_error,0.001); % define scale of bound circulation
    error = max(abs(GammaNew-Gamma)); % difference betweeen iterations
    error = error/ref_error; % relative error
    
    if error < errorlimit
        %kiter=Niterations;  % if error smaller than limit, stop iteration cycle
        disp(['Number of iterations to convergence: ',num2str(i)]);
        break
    end
    
    % set new estimate of bound circulation
    GammaNew = 0.7*Gamma + 0.3*GammaNew;
    %GammaNew(ig) = (1-ConvWeight)*Gamma(ig) + ConvWeight*GammaNew(ig);
end % end iteration loop

if i==Niterations
    disp("Solution not converged!");
end
end