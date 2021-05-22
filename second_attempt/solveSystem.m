function [a, aline, r_R, Fnorm, Ftan, GammaNew, Alpha, Inflow] = solveSystem(InfluenceMatrix, RotorWakeSystem, Radius, Omega, wind)
    
    UMat = InfluenceMatrix.U;
    VMat = InfluenceMatrix.V;
    WMat = InfluenceMatrix.W;
    
    Uinf = wind(1);
    Vinf = wind(2);
    Winf = wind(3);
    
    NBlades = RotorWakeSystem.NBlades;
%     cp = RotorWakeSystem.cp;%for solving using no-normal flow condition
    bound = RotorWakeSystem.bound;
    Ncp = RotorWakeSystem.NpanelsPerBlade*NBlades;  % as many points as there are panels
    Nrings = RotorWakeSystem.NpanelsPerBlade*NBlades;   % as many horseshoes as there are panels
    GammaNew = ones(Ncp,1); % initial guess, defined in control points
    a = zeros(Ncp,1);
    aline = zeros(Ncp,1);
    r_R = zeros(Ncp,1); % location of control point
    Fnorm = zeros(Ncp,1);
    Ftan = zeros(Ncp,1);
    Niter = 1000;   % maximum number of performed iterations
    errorlimit = 0.0001;    % error tolerance for the iterative process
    
    for i=1:Niter
        Gamma = GammaNew;   % update bound circulation
        
        % calculate velocity, circulation and loads at the controlpoints
        for icp= 1:Ncp  % loop over all control points
            localcp = bound.centrepoint(:,icp);
            % determine radial position of the controlpoint;
            radialpos = bound.radialcentrepoint(icp);%[m]
            u=0; % initialize velocity at control point
            v=0;
            w=0;
            % multiply icp line of Matrix with vector of circulation Gamma to calculate velocity at controlpoint
            for jring = 1:Nrings % loop over all vortex rings
                u = u + UMat(icp,jring)*Gamma(jring); % axial component of velocity
                v = v + VMat(icp,jring)*Gamma(jring); % y-component of velocity
                w = w + WMat(icp,jring)*Gamma(jring); % z-component of velocity
            end
            
            Urot = cross([-Omega,0,0],localcp);
            Urel = [Uinf+u+Urot(1), Vinf+v+Urot(2), Winf+w+Urot(3)];
            azimdir = cross([-1/radialpos, 0, 0]  , localcp); % rotational direction
            Utan = dot(azimdir, Urel); % azimuthal direction
            Uaxial =  dot([1, 0, 0] , Urel); % axial velocity

            [fnorm , ftan, gamma, alpha, inflow] = loadBladeElement(Utan, Uaxial, radialpos/Radius);

            % new point of new estimate of circulation for the blade section
            GammaNew(icp) = gamma;

            % update output vector
            a(icp) = (-(u + Urot(1))/wind(1));
            aline(icp) = Utan/(radialpos*Omega)-1;
            r_R(icp) = radialpos/Radius;    % location of control point
            Fnorm(icp) = fnorm ;
            Ftan(icp) = ftan ;
            Alpha(icp) = alpha;
            Inflow(icp) = inflow;
        end % end loop control points

        % check convergence of solution
        ref_error = max(abs(GammaNew));
        ref_error = max(ref_error,0.001); % define scale of bound circulation
        error = max(abs(GammaNew-Gamma)); % difference betweeen iterations
        error = error/ref_error; % relative error

        if error < errorlimit   % if error smaller than limit, stop iteration cycle
            disp(['Number of iterations to convergence: ',num2str(i)]);
            break
        end
        
        % set new estimate of bound circulation
        GammaNew = 0.7*Gamma + 0.3*GammaNew;
        disp(['Iter count:',num2str(i)]);
    end % end iteration loop

    if i==Niter
        disp("Solution not converged!");
    end
end