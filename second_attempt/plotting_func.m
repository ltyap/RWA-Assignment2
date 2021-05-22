function plotting_func(windvel,Radius, N, NBlades, Omega, a, aline, r_R_cp, ct,cp,cq, Gamma_temp, alpha, inflow)
colors = ['r'; 'b'; 'g'; 'm'; 'y']; % define color scheme, more colors need to be added for larger analysis
inflow = rad2deg(inflow); % change inflow angle from rad to degrees
Gamma_temp = Gamma_temp/(pi*norm(windvel)^2/(NBlades*Omega));

%% a and a'
figure()
hold on
plot(r_R_cp(1:N), a(1:N), '-', 'Color', colors(1), 'DisplayName',sprintf('a'))
plot(r_R_cp(1:N), aline(1:N), '--', 'Color', colors(2), 'DisplayName',sprintf('a'''))
title('Axial and tangential induction')
grid on
grid minor
xlabel('r/R')
legend show
hold off

%% alpha
figure()
hold on
plot(r_R_cp(1:N), alpha(1:N), '-', 'Color', colors(1))
title('Angle of attack \alpha [\circ]')
grid on
grid minor
xlabel('r/R')
% legend show
hold off

%% inflow
figure()
hold on
plot(r_R_cp(1:N), inflow(1:N), '-', 'Color', colors(1))
title('Angle of inflow [\circ]')
grid on
grid minor
xlabel('r/R')
% legend show
hold off

%% CT and CP
figure()
hold on
plot(r_R_cp(1:N), ct(1:N), '-', 'Color', colors(1), 'DisplayName',sprintf('C_T'))
plot(r_R_cp(1:N), cp(1:N), '--', 'Color', colors(2), 'DisplayName',sprintf('C_P'))
title('C_T and C_P')
grid on
grid minor
xlabel('r/R')
legend show
hold off

%% CQ
figure()
hold on
plot(r_R_cp(1:N), cq(1:N), '-', 'Color', colors(1))
title('C_Q')
grid on
grid minor
xlabel('r/R')
% legend show
hold off

%% CN

%% Circulation
figure()
hold on
plot(r_R_cp(1:N), Gamma_temp(1:N), '-', 'Color', colors(1))
title('Circulation distribution (non-dimensionalized by \pi U_\infty^2 / \Omega N_{blades})')
grid on
grid minor
xlabel('r/R')
% legend show
hold off
end