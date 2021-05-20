function [chord_distribution,twist_distribution] = BladeGeometry(r_R)
%Returns the chord and twist distribution along the blade

pitch = 2; % degrees
chord_distribution = 3*(1-r_R)+1; % [m]
twist_distribution = -14*(1-r_R)+pitch; % [deg]
end