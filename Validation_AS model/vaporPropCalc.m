function [rho_f, mu_f, k_f, cp_f] = vaporPropCalc(t_ref, y_ref, p)


R = 8.31447;
mw = 18.02e-3;
Rm = R/mw;
rho_f = (y_ref*p)/Rm/t_ref;

% mu_f = 10.6e-6;
Teske = [11.602, 11.600, 11.627, 11.603, 11.561,11.588, 11.596, 11.537, 11.527, 11.543];
mu_f = mean(Teske)*1e-6;
% mu_f = 0.5535e-3;

k_f = 23e-3;
% k_f = 22.7e-3;

% a = 32.24;
% b = 0.1923e-2;
% c = 1.055e-5;
% d = -3.595e-9;
% cp_f = a + b*t_ref + c*t_ref^2 + d*t_ref^3;
% cp_f = cp_f/18.02*1e3;
cp_f = 1.98e3;


end
















