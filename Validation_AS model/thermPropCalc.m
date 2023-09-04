function [pSat, l_Tsurf] = thermPropCalc(t_surf)

a = -7.77224;
b = 1.45684;
c = -2.71942;
d = -1.41336;
Pc = 221.2e5;
Tc = 647.3;
tau = 1 - t_surf/Tc;

pSat = exp((a*tau + b*tau^1.5 + c*tau^3 + d*tau^6)/(1 - tau))*Pc;

l_Tsurf = 2308e3;
% l_Tsurf = 1.91846e6*(t_surf/(t_surf - 33.91))^2;

end