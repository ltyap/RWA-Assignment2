function [r_R] = RadialSpacing(N, TipLocation_R, RootLocation_R, uniform)
%RADIALSPACING: Defines spacing for rotor blade
%Uniform Spacing
    if uniform
        r_R = linspace(RootLocation_R,TipLocation_R,N+1);
    else
        %cosine spacing
        theta = linspace(0,pi,N+1);
        b = (TipLocation_R-RootLocation_R)/2;
        r_R = RootLocation_R + flip(b*(cos(theta)+1)); % cosine spacing              
    end
end