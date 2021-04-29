% 3D velocity induced by a vortex filament
function u_ind = induced_v_from_vortex(gamma, point1, point2, cp)
% function to calculate the velocity induced by a straight 3D vortex filament
% with circulation GAMMA at a point VP1. The geometry of the vortex filament
% is defined by its edges: the filament starts at [X1, Y1, Z1] and ends at [X2, Y2, Z2].
% [XP, YP, ZP] denote the coordinates of target point where the velocity is calculated
% the input CORE defines a vortex core radius, inside which the velocity
% is defined as a solid body rotation.
% The function is adapted from the algorithm presented in:
%                Katz, Joseph, and Allen Plotkin. Low-speed aerodynamics.
%                Vol. 13. Cambridge university press, 2001.


% calculate geometric relations for integral of the velocity induced by filament
r1 = sqrt(sum((cp-point1).^2));
r2 = sqrt(sum((cp-point2).^2));
r12x = (cp(2)-point1(2))*(cp(3)-point2(3))-(cp(3)-point1(3))*(cp(2)-point2(2));
r12y = -(cp(1)-point1(1))*(cp(3)-point2(3))+(cp(3)-point1(3))*(cp(1)-point2(1));
r12z = (cp(1)-point1(1))*(cp(2)-point2(2))-(cp(2)-point1(2))*(cp(1)-point2(1));
r12_sqr = r12x^2+r12y^2+r12z^2;

r01 = sum((point2-point1).*(cp-point1));
r02 = sum((point2-point1).*(cp-point2));

% determine scalar
if r12_sqr==0
    K=0;
elseif r1==0
    K = gamma/(4*pi*r12_sqr)*(-r02/r2);
elseif r2==0
    K = gamma/(4*pi*r12_sqr)*(r01/r1);
else
    K = gamma/(4*pi*r12_sqr)*(r01/r1-r02/r2);   
end

% induced velocity - from bound vortex
u_ind = K*[r12x, r12y, r12z];
end