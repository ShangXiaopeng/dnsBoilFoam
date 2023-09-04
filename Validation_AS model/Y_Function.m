function dY = Y_Function(tt,Y,C_D,Rho_Inf,Rho_L,U_Inf,mDot,Q_L,Cp_L)

dY = zeros(3,1);
dY(1,1) = 3*C_D/8/Y(2)*Rho_Inf/Rho_L*abs(U_Inf - Y(1))*(U_Inf - Y(1));
dY(2,1) = -mDot/4/pi/Rho_L/Y(2).^2;
dY(3,1) = 3*Q_L/4/pi/Y(2).^3/Rho_L/Cp_L;
end