function [r_R, theta, b] = RadialSpacing(N, TipLocation_R, RootLocation_R)
%RADIALSPACING: Defines spacing for rotor blade
% cosine spacing for now

%cosine spacing
theta = linspace(0,pi,N+1);
b = (TipLocation_R-RootLocation_R)/2;
r_R = RootLocation_R + flip(b*(cos(theta)+1)); % cosine spacing
end