function [result] = Rotate(input,angle)

%rotation about x-axis;

A = [1,0,0;...
    0,cos(angle),-sin(angle);...
    0,sin(angle),cos(angle)];

result = A*input;
end

