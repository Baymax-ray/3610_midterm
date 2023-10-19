%% parameters
P_tot=107000;
V=3500;

DinM=30;% 30 days in a month
iniP_G1=P_tot*0.137; %<18
iniP_G2=P_tot*0.164; %>65
iniP_G3=P_tot-iniP_G1-iniP_G2; %18-65
%death rate
DR_G1=0.001;
DR_G2=0.18;
DR_G3=0.02;
%infection rate
IR_G1toG1=0.18*2; 
IR_G1toG2=0.28*2; 
IR_G1toG3=0.18*2; 
IR_G2toG1=0.17*2; 
IR_G2toG2=0.25*2; 
IR_G2toG3=0.17*2; 
IR_G3toG1=0.08*2; %the ratio of parents to children : children to parents is about 3:7
IR_G3toG2=0.28*2; 
IR_G3toG3=0.12*2; 


%recover days
R_G=3; %assume everyone is no longer contagious after 3 days
%VDR---vaccine effect on death rate
%VIR---vaccine effect on infection rate
VDR=0.1;
VIR=0.3;
FLW = 436; %unvaccinated front line workers
%% naive
[I_G1,I_G2,I_G3,V_I_G1,V_I_G2,V_I_G3,D_G1,D_G2,D_G3,V_D_G1,V_D_G2,V_D_G3,NV_P_G1,NV_P_G2,NV_P_G3]=naive(VDR,VIR,V,DinM,iniP_G1,iniP_G2,iniP_G3,DR_G1,DR_G2,DR_G3,IR_G1toG1,IR_G1toG2,IR_G1toG3,IR_G2toG1,IR_G2toG2,IR_G2toG3,IR_G3toG1,IR_G3toG2,IR_G3toG3,R_G);

%% front line workers+123
%[I_G1,I_G2,I_G3,V_I_G1,V_I_G2,V_I_G3,D_G1,D_G2,D_G3,V_D_G1,V_D_G2,V_D_G3,NV_P_G1,NV_P_G2,NV_P_G3]=seq123(VDR,VIR,V,DinM,iniP_G1,iniP_G2,iniP_G3,DR_G1,DR_G2,DR_G3,IR_G1toG1,IR_G1toG2,IR_G1toG3,IR_G2toG1,IR_G2toG2,IR_G2toG3,IR_G3toG1,IR_G3toG2,IR_G3toG3,R_G,FLW);
%% front line workers+213
%[I_G1,I_G2,I_G3,V_I_G1,V_I_G2,V_I_G3,D_G1,D_G2,D_G3,V_D_G1,V_D_G2,V_D_G3,NV_P_G1,NV_P_G2,NV_P_G3]=seq213(VDR,VIR,V,DinM,iniP_G1,iniP_G2,iniP_G3,DR_G1,DR_G2,DR_G3,IR_G1toG1,IR_G1toG2,IR_G1toG3,IR_G2toG1,IR_G2toG2,IR_G2toG3,IR_G3toG1,IR_G3toG2,IR_G3toG3,R_G,FLW);
%% ratio12
ra12=1000;
%[I_G1,I_G2,I_G3,V_I_G1,V_I_G2,V_I_G3,D_G1,D_G2,D_G3,V_D_G1,V_D_G2,V_D_G3,NV_P_G1,NV_P_G2,NV_P_G3]=ratio12(VDR,VIR,V,DinM,iniP_G1,iniP_G2,iniP_G3,DR_G1,DR_G2,DR_G3,IR_G1toG1,IR_G1toG2,IR_G1toG3,IR_G2toG1,IR_G2toG2,IR_G2toG3,IR_G3toG1,IR_G3toG2,IR_G3toG3,R_G,FLW,ra12);
%% 
I_G1=I_G1(1:50);
I_G2=I_G2(1:50);
I_G3=I_G3(1:50);
V_I_G1=V_I_G1(1:50);
V_I_G2=V_I_G2(1:50);
V_I_G3=V_I_G3(1:50);
D_G1=D_G1(1:50);
D_G2=D_G2(1:50);
D_G3=D_G3(1:50);
V_D_G1=V_D_G1(1:50);
V_D_G2=V_D_G2(1:50);
V_D_G3=V_D_G3(1:50);
NV_P_G1=NV_P_G1(1:50);
NV_P_G2=NV_P_G2(1:50);
NV_P_G3=NV_P_G3(1:50);


%% plot result
% plot infected population per day
figure(1)
plot(I_G1+I_G2+I_G3+V_I_G1+V_I_G2+V_I_G3,"Color",[83 81 84]./255,'LineWidth',2)
hold on
plot(I_G1+V_I_G1,"Color",[57 106 177]./255,'LineWidth',2)
plot(I_G2+V_I_G2,"Color",[204 37 41]./255,'LineWidth',2)
plot(I_G3+V_I_G3,"Color",[62 150 81]./255,'LineWidth',2)
xlabel('days')
ylabel('population')
title('new infected population at each day')
legend('total infected population','infected population in group 1','infected population in group 2','infected population in group 3')
hold off

% plot dead population per day
figure(2)
plot(D_G1+D_G2+D_G3+V_D_G1+V_D_G2+V_D_G3,"Color",[83 81 84]./255,"LineWidth",2)
hold on
plot(D_G1+V_D_G1,"Color",[57 106 177]./255,"LineWidth",2)
plot(D_G2+V_D_G2,"Color",[204 37 41]./255,"LineWidth",2)
plot(D_G3+V_D_G3,"Color",[62 150 81]./255,"LineWidth",2)
xlabel('days')
ylabel('population')
title('new dead population at each day')
legend('total dead population','dead population in group 1','dead population in group 2','dead population in group 3')
hold off

% plot accumulated infected population
figure(3)
plot(cumsum(I_G1+I_G2+I_G3+V_I_G1+V_I_G2+V_I_G3),"Color",[83 81 84]./255,"LineWidth",2)
hold on
plot(cumsum(I_G1+V_I_G1),"Color",[57 106 177]./255,"LineWidth",2)
plot(cumsum(I_G2+V_I_G2),"Color",[204 37 41]./255,"LineWidth",2)
plot(cumsum(I_G3+V_I_G3),"Color",[62 150 81]./255,"LineWidth",2)
xlabel('days')
ylabel('population')
title('accumulated infected population')
legend('total infected population','infected population in group 1','infected population in group 2','infected population in group 3')
hold off

% plot accumulated dead population
figure(4)
plot(cumsum(D_G1+D_G2+D_G3+V_D_G1+V_D_G2+V_D_G3),"Color",[83 81 84]./255,"LineWidth",2)
hold on
plot(cumsum(D_G1+V_D_G1),"Color",[57 106 177]./255,"LineWidth",2)
plot(cumsum(D_G2+V_D_G2),"Color",[204 37 41]./255,"LineWidth",2)
plot(cumsum(D_G3+V_D_G3),"Color",[62 150 81]./255,"LineWidth",2)
xlabel('days')
ylabel('population')
title('accumulated dead population')
legend('total dead population','dead population in group 1','dead population in group 2','dead population in group 3')
hold off

% plot vaccinated population per day
figure(5)
plot(NV_P_G1+NV_P_G2+NV_P_G3,"Color",[83 81 84]./255,"LineWidth",2)
hold on
plot(NV_P_G1,"Color",[57 106 177]./255,"LineWidth",2)
plot(NV_P_G2,"Color",[204 37 41]./255,"LineWidth",2)
plot(NV_P_G3,"Color",[62 150 81]./255,"LineWidth",2)
xlabel('days')
ylabel('population')
title('new vaccinated population')
legend('total vaccinated population', 'vaccinated population in group 1', 'vaccinated population in group2', 'vaccinated population in group3')
hold off

% plot accumulated vaccinated population
figure(6)
plot(cumsum(NV_P_G1+NV_P_G2+NV_P_G3),"Color",[83 81 84]./255,"LineWidth",2)
hold on
plot(cumsum(NV_P_G1),"Color",[57 106 177]./255,"LineWidth",2)
plot(cumsum(NV_P_G2),"Color",[204 37 41]./255,"LineWidth",2)
plot(cumsum(NV_P_G3),"Color",[62 150 81]./255,"LineWidth",2)
xlabel('days')
ylabel('population')
title('accumulated vaccinated population')
legend('total vaccinated population', 'vaccinated population in group 1', 'vaccinated population in group2', 'vaccinated population in group3')
hold off

%% 
t1=cumsum(D_G1+D_G2+D_G3+V_D_G1+V_D_G2+V_D_G3);
disp(t1(end));
t2=cumsum(I_G1+I_G2+I_G3+V_I_G1+V_I_G2+V_I_G3);
disp(t2(end));
total_death_G1=cumsum(D_G1+V_D_G1);
total_death_G2=cumsum(D_G2+V_D_G2);
total_death_G3=cumsum(D_G3+V_D_G3);
LOL=total_death_G1(end)*expect_lossG1+total_death_G2(end)*expect_lossG2+total_death_G3(end)*expect_lossG3