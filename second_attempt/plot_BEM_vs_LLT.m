%% plotting BEM vs LLT
clc
clear all
close all

windvel = [10,0,0];
TSR=8;
Radius=50;
Omega=norm(windvel)*TSR/Radius;
NBlades=3;

load('results_llt_N20.mat')
llt_N = length(results_llt)/NBlades;
load('results_bem_N20.mat');
bem_N = length(results_bem);

colors = ['r'; 'b'; 'g'; 'm'; 'y']; % define color scheme, more colors need to be added for larger analysis
% change inflow angle from rad to degrees
results_bem(:,9) = rad2deg(results_bem(:,9));
results_llt(:,9) = rad2deg(results_llt(:,9));
% non-dimensionalize circulation - only for llt!
results_llt(:,7) = results_llt(:,7)/(pi*norm(windvel)^2/(NBlades*Omega));

%% a
figure()
hold on
plot(results_llt(1:llt_N,1), results_llt(1:llt_N,2), '-', 'Color', colors(1), 'DisplayName',sprintf('LLT'))
plot(results_bem(:,1), results_bem(:,2), '--', 'Color', colors(2), 'DisplayName',sprintf('BEM'))
title('Axial induction')
grid on
grid minor
xlabel('r/R')
legend show
hold off

%% a'
figure()
hold on
plot(results_llt(1:llt_N,1), results_llt(1:llt_N,3), '-', 'Color', colors(1), 'DisplayName',sprintf('LLT'))
plot(results_bem(:,1), results_bem(:,3), '--', 'Color', colors(2), 'DisplayName',sprintf('BEM'))
title('Tangential induction')
grid on
grid minor
xlabel('r/R')
legend show
hold off

%% alpha
figure()
hold on
plot(results_llt(1:llt_N,1), results_llt(1:llt_N,8), '-', 'Color', colors(1), 'DisplayName',sprintf('LLT'))
plot(results_bem(:,1), results_bem(:,8), '--', 'Color', colors(2), 'DisplayName',sprintf('BEM'))
title('Angle of attack \alpha [\circ]')
grid on
grid minor
xlabel('r/R')
legend show
hold off

%% inflow
figure()
hold on
plot(results_llt(1:llt_N,1), results_llt(1:llt_N,9), '-', 'Color', colors(1), 'DisplayName',sprintf('LLT'))
plot(results_bem(:,1), results_bem(:,9), '--', 'Color', colors(2), 'DisplayName',sprintf('BEM'))
title('Angle of inflow [\circ]')
grid on
grid minor
xlabel('r/R')
legend show
hold off

%% CT
figure()
hold on
plot(results_llt(1:llt_N,1), NBlades*results_llt(1:llt_N,4), '-', 'Color', colors(1),'DisplayName',sprintf('LLT'))
plot(results_bem(:,1), results_bem(:,4), '--', 'Color', colors(2), 'DisplayName',sprintf('BEM'))
title('C_T')
grid on
grid minor
xlabel('r/R')
legend show
hold off

%% CP
figure()
hold on
plot(results_llt(1:llt_N,1), NBlades*results_llt(1:llt_N,5), '-', 'Color', colors(1), 'DisplayName',sprintf('LLT'))
plot(results_bem(:,1), results_bem(:,5), '--', 'Color', colors(2), 'DisplayName',sprintf('BEM'))
title('C_P')
grid on
grid minor
xlabel('r/R')
legend show
hold off

%% CQ
figure()
hold on
plot(results_llt(1:llt_N,1), NBlades*results_llt(1:llt_N,6), '-', 'Color', colors(1), 'DisplayName',sprintf('LLT'))
plot(results_bem(:,1), results_bem(:,6), '--', 'Color', colors(2), 'DisplayName',sprintf('BEM'))
title('C_Q')
grid on
grid minor
xlabel('r/R')
legend show
hold off

%% CN

%% Circulation
figure()
hold on
plot(results_llt(1:llt_N,1), results_llt(1:llt_N,7), '-', 'Color', colors(1), 'DisplayName',sprintf('LLT'))
plot(results_bem(:,1), results_bem(:,7), '--', 'Color', colors(2), 'DisplayName',sprintf('BEM'))
title('Circulation distribution (non-dimensionalized by \pi U_\infty^2 / \Omega N_{blades})')
grid on
grid minor
xlabel('r/R')
legend show
hold off

%% Fnorm
figure()
hold on
plot(results_llt(1:llt_N,1), results_llt(1:llt_N,10), '-', 'Color', colors(1), 'DisplayName',sprintf('LLT'))
plot(results_bem(:,1), results_bem(:,10), '--', 'Color', colors(2), 'DisplayName',sprintf('BEM'))
title('Fnorm')
grid on
grid minor
xlabel('r/R')
legend show
hold off

%% Ftan
figure()
hold on
plot(results_llt(1:llt_N,1), results_llt(1:llt_N,11), '-', 'Color', colors(1), 'DisplayName',sprintf('LLT'))
plot(results_bem(:,1), results_bem(:,11), '--', 'Color', colors(2), 'DisplayName',sprintf('BEM'))
title('Ftan')
grid on
grid minor
xlabel('r/R')
legend show
hold off