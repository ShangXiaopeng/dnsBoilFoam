clear all
clc

% % % % % % % Appendix B _ Wiley % % % % % % % %

Rho_L = 971.8;
Cp_L = 4194;
K_L = 0.6562;
Mu_L = 0.0003545;

Alpha_L = K_L/Rho_L/Cp_L;

Sigma_L = 0.06267;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ar =1/3;

Y_Inf = 0;

MW_F = 18.02e-3;
MW_A = 28.97e-3;

P = 101325;

T_Inf = 363;

Rho_Inf = rhoInfCalc(T_Inf, P);

U_Inf = 0.0;

RT = 2;

dt = 1e-3;
t = 0:dt:RT;
N = length(t);

U = zeros(1,N);
R_Surf = zeros(1,N);
T_Surf = zeros(1,N);
mDot = zeros(1,N);
x = zeros(1,N);
x(1) = 0;
U(1) = 10e-3;
R_Surf(1) = 1.0e-3;
T_Surf(1) = 343;

Q_L = 0;
% B_T = 1;


% % % % % % % % % % % % % % % % % % % % % % % 

NR = 100;
T_Old = zeros(1,NR+1);
T_New = zeros(1,NR+1);
T_Old(1:NR+1) = T_Surf(1);

% % % % % % % % % % % % % % % % % % % % % % % % % % % 


for kk = 1:N

    epsilon_T = 0.01;
    epsilon_mDot = 0.01;
    
    absTol_T = inf;
    relTol_mDot = inf;
    
    U_Pre = U(kk);
    R_Surf_Pre = R_Surf(kk);
    T_Surf_Pre = T_Surf(kk);
    
    count = 1;
    
    while (absTol_T > epsilon_T)||(relTol_mDot > epsilon_mDot)
        
        if count <= 100
            
            [pSat_Surf, L_Surf] = thermPropCalc(T_Surf_Pre);
            
            Y_Surf = (pSat_Surf/P)*MW_F/(pSat_Surf/P*MW_F + (1-pSat_Surf/P)*MW_A);
        
            T_Ref = T_Surf_Pre + Ar*(T_Inf - T_Surf_Pre);
        
            Y_Ref = Y_Surf + Ar*(Y_Inf - Y_Surf);
        
            [Rho_F, Mu_F, K_F, Cp_F] = vaporPropCalc(T_Ref, Y_Ref, P);
        
            [Rho_G, Cp_G, Mu_G, K_G, D_G, Le_G, Pr_G, Sc_G] = mixturePropCalc(T_Ref, P, Y_Ref);
        
            Re = 2*Rho_Inf*abs(U_Pre - U_Inf)*R_Surf_Pre/Mu_G;
        
            Nu_0 = 2 + 0.552*Re^(1/2)*Pr_G^(1/3);
            Sh_0 = 2 + 0.552*Re^(1/2)*Sc_G^(1/3);
        
%             Nu_0 = 1 + (1 + Re*Pr_G)^(1/3)*(1*(Re<=1) + Re^(0.077)*(Re>1));      
%             Sh_0 = 1 + (1 + Re*Sc_G)^(1/3)*(1*(Re<=1) + Re^(0.077)*(Re>1));
        
            B_M = (Y_Surf - Y_Inf)/(1 - Y_Surf); 
            F_M = (1 + B_M)^(0.7)*log(1 + B_M)/B_M;
            
            Sh_Asterisk = 2 + (Sh_0 - 2)/F_M;
        
            if count == 1
                
                mDot(kk) = 2*pi*Rho_G*D_G*R_Surf(kk)*Sh_Asterisk*log(1 + B_M);
                mDot_Pre = mDot(kk);
                
            else
                
                mDot_Pre = 2*pi*Rho_G*D_G*R_Surf_Pre*Sh_Asterisk*log(1 + B_M);
                
            end
            
            B_T = Cp_F*(T_Inf - T_Surf_Pre)/(L_Surf + Q_L/mDot_Pre);

            absTol_B = inf;
            epsilon_B = 1e-5;

            while (absTol_B > epsilon_B)
            
                F_T = (1 + B_T)^(0.7)*log(1 + B_T)/B_T;           
                Nu_Asterisk = 2 + (Nu_0 - 2)/F_T;
 
                Phi = (Cp_F/Cp_G)*(Sh_Asterisk/Nu_Asterisk)*(1/Le_G);
                B_T_New = (1 + B_M)^(Phi) - 1;
            
                absTol_B = abs(B_T - B_T_New);
                B_T = B_T_New;
        
            end
            
            Q_L = mDot_Pre*(Cp_F*(T_Inf - T_Surf_Pre)/B_T - L_Surf);
            C_D = 24/Re*(1 + Re^(2/3)/6);


% % % % % % % % % % % % % % % % % % update U, Rs and mdot % % % % % % % % % % % % % % % % % % % %
            U_New = U(kk) + dt*(3*C_D/2/R_Surf(kk)*Rho_Inf/Rho_L*abs(U_Inf - U(kk))*(U_Inf - U(kk)));
            R_Surf_New = R_Surf(kk) + dt*(-mDot(kk)/4/pi/Rho_L/R_Surf(kk)^2);
            mDot_New = 2*pi*Rho_G*D_G*R_Surf(kk)*Sh_Asterisk*log(1 + B_M);


%             U_New = U(kk) + dt*(3*C_D/2/R_Surf_Pre*Rho_Inf/Rho_L*abs(U_Inf - U_Pre)*(U_Inf - U_Pre));
%             R_Surf_New = R_Surf(kk) + dt*(-mDot_Pre/4/pi/Rho_L/R_Surf_Pre^2);
%             mDot_New = 2*pi*Rho_G*D_G*R_Surf_Pre*Sh_Asterisk*log(1 + B_M);

% % % % % % % % % % % % % % % % % % update U, Rs and mdot % % % % % % % % % % % % % % % % % % % %


% % % % % % % % % % % % % % % % % Rapid Mixing Limit % % % % % % % % % % % % % % % % % % % % 

            T_Surf_New = T_Surf(kk) + dt*(3*Q_L/4/pi/R_Surf(kk)^3/Rho_L/Cp_L);

            
%             T_Surf_New = T_Surf(kk) + dt*(3*Q_L/4/pi/R_Surf_Pre^3/Rho_L/Cp_L); % not recommended  

% % % % % % % % % % % % % % % % % Rapid Mixing Limit % % % % % % % % % % % % % % % % % % % % 



% % % % % % % % % % % % % % % % Conduction Limit % % % % % % % % % % % % % % % % % 
% 
%             deltaR = R_Surf_Pre/NR;
%             
%             for i1 = 2:NR
%                 
%                 r = 0 + (i1-1)*deltaR;
% 
%                 AA(i1) = 1 + 2*Alpha_L*dt/deltaR^2;
%                 BB(i1) = Alpha_L*dt*(1/deltaR^2 + 1/r/deltaR);
%                 CC(i1) = Alpha_L*dt*(1/deltaR^2 - 1/r/deltaR);
%                 DD(i1) = T_Old(i1);
%                 
%             end
% 
%             AA(1) = 1;
%             BB(1) = 1;
%             CC(1) = 0;
%             DD(1) = 0;
%             
%             AA(NR+1) = 1;
%             BB(NR+1) = 0;
%             CC(NR+1) = 1;
%             DD(NR+1) = Q_L*deltaR/4/pi/K_L/R_Surf_Pre^2; 
%             
%             T_New = TDMA(AA, BB, CC, DD, NR+1);
%             T_Surf_New = T_New(NR+1);
% 
% % % % % % % % % % % % % % % % % Conduction Limit % % % % % % % % % % % % % % % % % 




% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
            
            absTol_T = abs(T_Surf_New - T_Surf_Pre);
            relTol_mDot = abs((mDot_New - mDot_Pre)/mDot_Pre);
                
            U_Pre = (U_New + U_Pre)/2;
            R_Surf_Pre = (R_Surf_New + R_Surf_Pre)/2;
            T_Surf_Pre = (T_Surf_New + T_Surf_Pre)/2;
            mDot_Pre = (mDot_New + mDot_Pre)/2;
            
       
        else
            
            disp('Unable to get converged in 100 iterations!');
            
        end
        
        count = count + 1;
        
    end

    U(kk+1) = U_New;
    R_Surf(kk+1) = R_Surf_New;
    T_Surf(kk+1) = T_Surf_New;
    mDot(kk+1) = mDot_New;

    x(kk+1) = x(kk) + U(kk)*dt;
    disp(t(kk));

    T_Old = T_New;
end



plot(t(1:kk),R_Surf(1:kk).^2/R_Surf(1)^2,'LineWidth',1.0,'MarkerSize',6.0)

set(gca,'FontSize',12,'TickLength',[0.02 0.02])


tt = [
    0.0000
    0.1000
    0.2000
    0.3000
    0.4000
    0.5000
    0.6000
    0.7000
    0.8000
    0.9000
    1.0000
    1.1000
    1.2000
    1.3000
    1.4000
    1.5000
    1.6000
    1.7000
    1.8000
    1.9000    
];

VW = 1e-8*[
   0.423896000000000
   0.422909000000000
   0.422073000000000
   0.421234000000000
   0.420461000000000
   0.419649000000000
   0.418897000000000
   0.418217000000000
   0.417519000000000
   0.416917000000000
   0.416205000000000
   0.415574000000000
   0.414967000000000
   0.414277000000000
   0.413711000000000
   0.413159000000000
   0.412627000000000
   0.412099000000000
   0.411576000000000
   0.411009000000000
];


RW = (VW*3/4/pi).^(1/3);

hold on 

plot(tt, RW.^2/RW(1)^2,'s','LineWidth',1.0,'MarkerSize',7.0)


xtick = 0.0:0.5:2.0;
set(gca,'XTick',xtick,'XTickLabel',{'0.0','0.5', '1.0', '1.5', '2.0'},'FontName','Times New Roman')
ytick = 0.97:0.01:1.0;
set(gca,'YTick',ytick,'YTickLabel',{'0.97', '0.98', '0.99',  '1.00'},'FontName','Times New Roman')

xlabel('\fontsize{11}TIME (ms)','FontName','Times New Roman')
% ylabel('\fontsize{12}NONDIMENSIONAL RADIUS')
ylabel('\fontsize{12}R^{2}/R^{2}_{0}','FontName','Times New Roman')

axis([0 2 0.97 1.0001])
% legend('Rapid Mixing Limit', 'Conduction Limit')

set(gcf,'Units','centimeters','Position',[1 2 17.5 15]);
set(gca,'Position',[0.175 0.17 0.775 0.78])


