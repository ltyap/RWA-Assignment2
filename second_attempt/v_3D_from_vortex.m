% 3D velocity induced by a vortex filament
function u_ind = v_3D_from_vortex(gamma, point1, point2, cp, CORE)
% calculate geometric relations for integral of the velocity induced by filament
point1 = point1(:);
point2 = point2(:);
cp = cp(:);
r1 = sqrt(sum((cp-point1).^2));
r2 = sqrt(sum((cp-point2).^2));

X1 = point1(1);
Y1 = point1(2);
Z1 = point1(3);

X2 = point2(1);
Y2 = point2(2);
Z2 = point2(3);

XP = cp(1);
YP = cp(2);
ZP = cp(3);

r1xr2_x = (YP-Y1)*(ZP-Z2)-(ZP-Z1)*(YP-Y2);
r1xr2_y = -(XP-X1)*(ZP-Z2)+(ZP-Z1)*(XP-X2);
r1xr2_z = (XP-X1)*(YP-Y2)-(YP-Y1)*(XP-X2);
r1xr_sqr = r1xr2_x^2+r1xr2_y^2+r1xr2_z^2;

r0r1 = (X2-X1)*(XP-X1)+(Y2-Y1)*(YP-Y1)+(Z2-Z1)*(ZP-Z1);
r0r2 = (X2-X1)*(XP-X2)+(Y2-Y1)*(YP-Y2)+(Z2-Z1)*(ZP-Z2);

% determine scalar
if (r1xr_sqr < CORE^2)
    r1xr_sqr=CORE^2;
end

if (r1 < CORE)
    r1 = CORE;
end

if (r2 < CORE)
    r2 = CORE;
end
%scalar
K = gamma/(4*pi*r1xr_sqr)*(r0r1/r1-r0r2/r2); 
% induced velocity
u_ind = K*[r1xr2_x, r1xr2_y, r1xr2_z];
u_ind = u_ind(:);
end