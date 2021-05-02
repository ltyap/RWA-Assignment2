% new iterative loop
function [a, aline, r_R, Fnorm, Ftan, Gamma_temp] = solveGamma(Uinf, N, Radius, cp, Influence_u, Influence_v, Influence_w, Omega, polar_alpha, polar_cl, polar_cd, cp_twist)
GammaNew = ones(N,1); % initial guess
U_inf = [Uinf,0,0];

Niterations = 1200;
errorlimit = 0.01;
error = 1.0;

% Niterations = 10000;    % maximum number of performed iterations
% Erroriterations = 0.0000001; % error tolerance for the iterative process

for i=1:Niterations
    Gamma = GammaNew;%???? % update bound circulation
    
    % calculate velocity, circulation and loads at the controlpoints
    for icp= 1:length(cp(:,1))
        
        % determine radial position of the controlpoint;
        r = cp(icp,2);  %   ??? not sure, mine
        %cp(icp,1) = 0;   % calculate everything wrt to c/4
        %r = sqrt(dot(cp(icp,:),cp(icp,:))); % this is from the tutorial
        
        u=0; % initialize velocity
        v=0;
        w=0;
        
        % multiply icp line of Matrix with vector of circulation Gamma to calculate velocity at controlpoint
        for j = 1:length(N)
            u = u + Influence_u(icp,j)*Gamma(j); % axial component of velocity
            v = v + Influence_v(icp,j)*Gamma(j); % y-component of velocity
            w = w + Influence_w(icp,j)*Gamma(j); % z-component of velocity
        end
        
        % calculate total perceived velocity
%         vrot = cross([-Omega, 0 , 0]  , cp(icp,:) ); % rotational velocity
%         vel1 = [U_inf(1)+ u + vrot(1), U_inf(2)+ v + vrot(2) , U_inf(3)+ w + vrot(3)]; % total perceived velocity at section
        
        % calculate azimuthal and axial velocity
%         azimdir = cross([-1/r, 0 , 0]  , cp(icp,:)); % rotational direction
%         Utan = dot(azimdir, vel1); % azimuthal direction
%         Uaxial =  dot([1, 0, 0] , vel1); % axial velocity
        
        %gives different results - not sure, mine
        Uaxial = Uinf+w;
        Utan = Omega*r+u;%dot([Uinf+u, v, w],[1,0,0]);
        Uper = sqrt(Uaxial^2+Utan^2); % for checking
        
        [fnorm , ftan, gamma, ~, ~] = loadBladeElement(r/Radius, cp(icp,1)*4, cp_twist(icp), polar_alpha, polar_cl, polar_cd, Uaxial, Utan);
        
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