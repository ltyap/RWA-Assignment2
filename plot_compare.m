function plot_compare(filenames, name_array, var_change, N, rotor, rho, windvel)
spacing_names={'cosine', 'uniform'}
NBlades=rotor.NBlades;
TSR=rotor.TSR;
Radius=rotor.Radius;
Omega = norm(windvel)*TSR/Radius;

for i=1:length(name_array)
    filename = filenames{i};
    temp = load(filename);
    results(:,:,i) = temp.results;
end

%% plotting
colors = ['r'; 'b'; 'g'; 'm'; 'k']; % define color scheme, more colors need to be added for larger analysis
styles = {'-', '--', '-.',':', '-','--'};
% change inflow angle from rad to degrees
results(:,9,:) = rad2deg(results(:,9,:));
% non-dimensionalize circulation - only for llt!
results(:,7,:) = results(:,7,:)/(pi*norm(windvel)^2/(NBlades*Omega));
% non-dimensionalize Fnorm and Ftan
results(:,10,:) = results(:,10,:)/(0.5*rho*norm(windvel)^2*Radius);
results(:,11,:) = results(:,11,:)/(0.5*rho*norm(windvel)^2*Radius);

%% a
figure()
hold on
for i=1:length(name_array)
    plot(results(1:N,1,i), results(1:N,2,i), 'linestyle',styles{i}, 'Color', colors(i), 'DisplayName',sprintf('%s', spacing_names{i}))
end
title('Axial induction')
grid on
grid minor
xlabel('r/R')
legend show
hold off

%% a'
figure()
hold on
for i=1:length(name_array)
    plot(results(1:N,1,i), results(1:N,3,i), 'linestyle',styles{i}, 'Color', colors(i), 'DisplayName',sprintf('%s', spacing_names{i}))
end
title('Tangential induction')
grid on
grid minor
xlabel('r/R')
legend show
hold off

%% alpha
figure()
hold on
for i=1:length(name_array)
    plot(results(1:N,1,i), results(1:N,8,i), 'linestyle',styles{i}, 'Color', colors(i), 'DisplayName',sprintf('%s', spacing_names{i}))
end
title('Angle of attack \alpha [\circ]')
grid on
grid minor
xlabel('r/R')
legend show
hold off

%% inflow
figure()
hold on
for i=1:length(name_array)
    plot(results(1:N,1,i), results(1:N,9,i), 'linestyle',styles{i}, 'Color', colors(i), 'DisplayName',sprintf('%s', spacing_names{i}))
end
title('Angle of inflow [\circ]')
grid on
grid minor
xlabel('r/R')
legend show
hold off

%% CT
figure()
hold on
for i=1:length(name_array)
    plot(results(1:N,1,i), results(1:N,4,i), 'linestyle',styles{i}, 'Color', colors(i), 'DisplayName',sprintf('%s', spacing_names{i}))
end
title('C_T')
grid on
grid minor
xlabel('r/R')
legend show
hold off

%% CP
figure()
hold on
for i=1:length(name_array)
    plot(results(1:N,1,i), results(1:N,5,i), 'linestyle',styles{i}, 'Color', colors(i), 'DisplayName',sprintf('%s', spacing_names{i}))
end
title('C_P')
grid on
grid minor
xlabel('r/R')
legend show
hold off

%% CQ
figure()
hold on
for i=1:length(name_array)
    plot(results(1:N,1,i), results(1:N,6,i), 'linestyle',styles{i}, 'Color', colors(i), 'DisplayName',sprintf('%s', spacing_names{i}))
end
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
for i=1:length(name_array)
    plot(results(1:N,1,i), results(1:N,7,i), 'linestyle',styles{i}, 'Color', colors(i), 'DisplayName',sprintf('%s', spacing_names{i}))
end
title('Circulation distribution (non-dimensionalized by \pi U_\infty^2 / \Omega N_{blades})')
grid on
grid minor
xlabel('r/R')
legend show
hold off

%% Fnorm
figure()
hold on
for i=1:length(name_array)
    plot(results(1:N,1,i), results(1:N,10,i), 'linestyle',styles{i}, 'Color', colors(i), 'DisplayName',sprintf('%s', spacing_names{i}))
end
title('Normal force (non-dimensionalized by 1/2 \rho U_\infty^2 R)')
grid on
grid minor
xlabel('r/R')
legend show
hold off

%% Ftan
figure()
hold on
for i=1:length(name_array)
    plot(results(1:N,1,i), results(1:N,11,i), 'linestyle',styles{i}, 'Color', colors(i), 'DisplayName',sprintf('%s', spacing_names{i}))
end
title('Tangential force (non-dimensionalized by 1/2 \rho U_\infty^2 R)')
grid on
grid minor
xlabel('r/R')
legend show
hold off
end