function [polar_alpha, polar_cl, polar_cd]=import_polars()
%% import polars
opts = delimitedTextImportOptions("NumVariables", 4);
opts.DataLines = [2, Inf];
opts.Delimiter = " ";
opts.VariableNames = ["Alpha", "Cl", "Cd", "Cm"];
opts.VariableTypes = ["double", "double", "double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";
opts.LeadingDelimitersRule = "ignore";
DU95W180 = readtable("DU95W180.txt", opts);
clear opts
polar_alpha = DU95W180.Alpha;
polar_cl = DU95W180.Cl;
polar_cd = DU95W180.Cd;

% %plot polars of the airfoil C-alfa and Cl-Cd
% figure()
% subplot(1,2,1);
% plot(polar_alpha,polar_cl);
% xlabel('\alpha');
% ylabel('C_{l}');
% xlim([-30 30]);
% grid on
% grid minor
% 
% subplot(1,2,2);
% plot(polar_alpha,polar_cd);
% xlabel('\alpha');
% ylabel('C_{d}');
% % xlim([0 0.1]);
% grid on
% grid minor
end