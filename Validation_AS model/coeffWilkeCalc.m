function phi = coeffWilkeCalc(mu, mw, ii, jj)

mui = mu(ii);
muj = mu(jj);

mwi = mw(ii);
mwj = mw(jj);

phi = (1 + (mui/muj)^0.5*(mwj/mwi)^0.25)^2/sqrt(8)/(1 + (mwi/mwj))^0.5;

end