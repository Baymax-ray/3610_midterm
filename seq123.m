function [I_G1,I_G2,I_G3,V_I_G1,V_I_G2,V_I_G3,D_G1,D_G2,D_G3,V_D_G1,V_D_G2,V_D_G3]=seq123(VDR,VIR,V,DinM,iniP_G1,iniP_G2,iniP_G3,DR_G1,DR_G2,DR_G3,IR_G1toG1,IR_G1toG2,IR_G1toG3,IR_G2toG1,IR_G2toG2,IR_G2toG3,IR_G3toG1,IR_G3toG2,IR_G3toG3,R_G)
% people always acccept vaccine
% prioritize kids>old>mid

%% initialize parameters
% unvaccinated uninfected population
P_G1=iniP_G1;
P_G2=iniP_G2;
P_G3=iniP_G3;
% vaccinated uninfected population
V_P_G1=0; %no one is vaccinated at the beginning
V_P_G2=0;
V_P_G3=0;
% unvaccinated infected population
I_G1=[5]; %initial infected population
I_G2=[5];
I_G3=[5];
% vaccinated infected population
V_I_G1=[0];
V_I_G2=[0];
V_I_G3=[0];
% unvaccinated dead population
D_G1=[0];
D_G2=[0];
D_G3=[0];
% vaccinated dead population
V_D_G1=[0];
V_D_G2=[0];
V_D_G3=[0];

end

