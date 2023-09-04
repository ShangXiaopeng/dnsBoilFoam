function rho_inf = rhoInfCalc(t_inf, p_inf)


R = 8.31447;
mw = 28.97e-3;
Rm = R/mw;
rho_inf = p_inf/Rm/t_inf;


end