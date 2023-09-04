function [rho_a, mu_a, k_a, cp_a] = airPropCalc(t_ref, y_ref, p)


R = 8.31447;
mw = 28.97e-3;
Rm = R/mw;
rho_a = (1 - y_ref)*p/Rm/t_ref;


mu0 = 1.716e-5;
mu_a = mu0*(273.11 + 110.56)/(t_ref + 110.56)*(t_ref/273.11)^1.5;

k_a = 1.5207e-11*t_ref^3 - 4.8574e-08*t_ref^2 + 1.0184e-4*t_ref - 3.9333e-4;
% k_a = 0.0302;
% cp_a = (0.251625 - 9.2525e-5*t_ref + 2.1334e-7*t_ref^2 - 1.0043e-10*t_ref^3)*4184;
cp_a = 1.0089e3;
end