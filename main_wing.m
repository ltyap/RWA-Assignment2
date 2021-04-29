% straight wing in steady flight with AoA=0;
clc
clear all
close all

Uinf=10;

% assume rectangular wing
N_segments=31;
chord=1;
AR = 66;
b = chord*AR;
S = b*chord;
c = chord*ones(1,N_segments);
alpha=deg2rad(5);  % in degrees
cl = 2*pi*sin(alpha);

% uniform distribution
% b_seg = b/N_segments*ones(1,N_segments); % segment lengths - assume

% cosine distribution
theta = linspace(0,pi,N_segments+1);
r = b/2;

% coordinates of the vortex filament endpoints
vfp(:,1) = [c/4, chord/4]; % x-coordinate
vfp(:,2) = flip(b/2*(cos(theta)+1)); % cosine spacing
% vfp(:2) = linspace(0,b,N_segments+1); % uniform spacing
vfp(:,3) = zeros(1,N_segments+1);

b_seg = vfp(2:end,2)-vfp(1:end-1,2);   % segment width
Gamma = 1/2*Uinf*cl*c';    % circulation at each segment

% inflow velocities at control points
vel = repmat(Uinf*[cos(alpha),0, sin(alpha)],N_segments,1);

% define control points at c/4
cp(:,1) = c/4; % x-coordinate
cp(:,3) = zeros(1,N_segments);     % z-coordinate

% % weird scaling from article - not middle of segment
% cp(1,2) = b_seg(1)^2/(b_seg(1)+b_seg(2))+vfp(1,2);
% cp(end,2) = b_seg(end-1)*b_seg(end)/(b_seg(end)+b_seg(end-1))+vfp(end-1,2);
% for i=2:length(b_seg)-1
%     cp(i,2) = 1/4*(b_seg(i-1)/(b_seg(i-1)+b_seg(i))+b_seg(i+1)/(b_seg(i+1)+b_seg(i))+1)*b_seg(i)+vfp(i,2);
% end

% at middle of each segment
cp(:,2) = vfp(1:end-1,2)+b_seg/2;     % y-coordinate


% vortex filaments at wing tips - assume length 1 convected in the x
% direction - THIS HAS TO BE ADDED INTO THE FUNCTION!
trail_vortex = 1e6; % trailing vortex length, almost infinity
% vfp_tip = [vfp(1,:); vfp(1,:);vfp(end,:);vfp(end,:)];
% vfp_tip(1,1) = vfp_tip(1,1)+trail_vortex;   % move the first one not second so signs agree
% vfp_tip(end,1) = vfp_tip(end,1)+trail_vortex;


% calculate induced velocity at each control point
% induced_vel = zeros(length(cp(:,1)),3);

for i=1:length(cp(:,1))
    % we know this will be zero for this case (rectangular wing), but not
    % for others
    for j=1:length(vfp(:,1))-1
        induced_vel=[0,0,0];
        temp = induced_v_from_vortex(1,vfp(j,:), vfp(j+1,:), cp(i,:));
        induced_vel=induced_vel+temp;
        
        % add influence of two trailing vortices
        temp = induced_v_from_vortex(1,vfp(j,:)+[trail_vortex,0,0], vfp(j,:), cp(i,:));
        induced_vel=induced_vel+temp;
        
        temp = induced_v_from_vortex(1,vfp(j+1,:), vfp(j+1,:)+[trail_vortex,0,0], cp(i,:));
        induced_vel=induced_vel+temp;
        
        % write into influence matrix - for unit gamma
        Influence_u(i,j) = induced_vel(1);
        Influence_v(i,j) = induced_vel(2);
        Influence_w(i,j) = induced_vel(3);
    end
end

% actual induced velocity (when we already have a proper gamma value)
actual_induced_vel(:,1) = Influence_u*Gamma;
actual_induced_vel(:,2) = Influence_v*Gamma;
actual_induced_vel(:,3) = Influence_w*Gamma;

% induced AoA
alpha_ind = atan(actual_induced_vel(:,3)/Uinf);
alpha = alpha+alpha_ind;

cl = 2*pi*sin(alpha);

plot(cp(:,2)/b, cl)
grid on