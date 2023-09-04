function [rho_g, cp_g, mu_g, k_g, d_g, Le_g, Pr_g, Sc_g] = mixturePropCalc(t_ref, p, y_f)


[rho_f, mu_f, k_f, cp_f] = vaporPropCalc(t_ref, y_f, p);
[rho_a, mu_a, k_a, cp_a] = airPropCalc(t_ref, y_f, p);


% rho_g = (y_f/rho_f + (1-y_f)/rho_a)^(-1);
% cp_g = cp_f *y_f + cp_a*(1-y_f);

mu = [mu_f, mu_a];
k = [k_f, k_a];
mw = [18.02e-3, 28.97e-3];
y = [y_f, 1-y_f];


mu_g = 0;
k_g = 0;
for i1 = 1:2
    
    sum_denominator = 0;
    for j1 = 1:2
        
        phi_ij = coeffWilkeCalc(mu, mw, i1, j1);
        sum_denominator = sum_denominator + y(j1)*phi_ij; 
    end

    mu_g = mu_g + y(i1)*mu(i1)/sum_denominator;
    k_g = k_g + y(i1)*k(i1)/sum_denominator;

end


rho_g = rho_f*y_f + rho_a*(1 - y_f);
cp_g = (rho_f*y_f*cp_f + rho_a*(1 - y_f)*cp_a)/rho_g;
% mu_g = (mu_f*y_f + mu_a*(1 - y_f));
% k_g = (k_f*y_f + k_a*(1 - y_f));








d_g = 3.5699e-5;
% d_g = 3.573e-5;
Le_g = k_g/rho_g/d_g/cp_g;
Pr_g = cp_g*mu_g/k_g;
Sc_g = mu_g/rho_g/d_g;
end





