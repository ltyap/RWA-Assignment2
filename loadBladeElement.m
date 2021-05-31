function [fnorm , ftan, gamma, alpha, inflowangle] = loadBladeElement(Utan, Uaxial, r_R)
% calculates the load in the blade element
vmag2 = Uaxial.^2 + Utan.^2;    % square of resultant velocity
inflowangle = atan2(Uaxial, Utan);  % calculate inflow angle
[chord, twist] = BladeGeometry(r_R);
alpha = twist + rad2deg(inflowangle);     % calculate angle of attack

[polar_alpha, polar_cl, polar_cd] = import_polars();

cl = interp1(polar_alpha, polar_cl, alpha,'linear','extrap');     % local cl from airfoil polars
cd = interp1(polar_alpha, polar_cd, alpha,'linear','extrap');     % local cd from airfoil polars

% if alpha  is out of the polar range - take maximum or minimum given value
% - not exactly physical
for i=1:length(alpha)
    if alpha(i) >= max(polar_alpha)
        [value, index] = max(polar_alpha);
        cl(i) = polar_cl(index);
        cd(i) = polar_cd(index);
    elseif alpha(i) <= min(polar_alpha)
        [value, index] = min(polar_alpha);
        cl(i) = polar_cl(index);
        cd(i) = polar_cd(index);
    end
end

lift = 0.5*vmag2.*cl*chord; % local lift
drag = 0.5*vmag2.*cd*chord; % local drag

fnorm = lift.*cos(inflowangle)+drag.*sin(inflowangle);  % local axial force
ftan = lift.*sin(inflowangle)-drag.*cos(inflowangle);   % local tangent force
gamma = 0.5*sqrt(vmag2).*cl*chord;  % local circulation
end