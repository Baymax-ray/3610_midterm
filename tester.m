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
VDR=0.999;
VIR=0.458;
FLW = 436; %unvaccinated front line workers


%potential_VIR=linspace(0.9,1,10);
potention_ra12=logspace(-4,4,40);
%for j=1:length(potential_VIR)
    %VIR=potential_VIR(j);
%%
total_death=zeros(length(potention_ra12),1);
total_death_G1=zeros(length(potention_ra12),1);
total_death_G2=zeros(length(potention_ra12),1);
total_death_G3=zeros(length(potention_ra12),1);
total_infected=zeros(length(potention_ra12),1);
for i=1:length(potention_ra12)
    
    ra12=potention_ra12(i);
    [I_G1,I_G2,I_G3,V_I_G1,V_I_G2,V_I_G3,D_G1,D_G2,D_G3,V_D_G1,V_D_G2,V_D_G3,NV_P_G1,NV_P_G2,NV_P_G3]=ratio12(VDR,VIR,V,DinM,iniP_G1,iniP_G2,iniP_G3,DR_G1,DR_G2,DR_G3,IR_G1toG1,IR_G1toG2,IR_G1toG3,IR_G2toG1,IR_G2toG2,IR_G2toG3,IR_G3toG1,IR_G3toG2,IR_G3toG3,R_G,FLW,ra12);
    cum_death=cumsum(V_D_G1+V_D_G2+V_D_G3+D_G1+D_G2+D_G3);
    cum_infected=cumsum(V_I_G1+V_I_G2+V_I_G3+I_G1+I_G2+I_G3);
    cum_death_G1=cumsum(V_D_G1+D_G1);
    cum_death_G2=cumsum(V_D_G2+D_G2);
    cum_death_G3=cumsum(V_D_G3+D_G3);
    total_death(i)=cum_death(end);
    total_death_G1(i)=cum_death_G1(end);
    total_death_G2(i)=cum_death_G2(end);
    total_death_G3(i)=cum_death_G3(end);
    total_infected(i)=cum_infected(end);
end

%% plot
% death and ratio
figure(1)
plot(potention_ra12,total_death,'LineWidth',2)
set(gca, 'XScale', 'log');
xlabel('ratio of kid:old')
ylabel('total death at the end')
title('total death vs ratio of kid:old')

% infected and ratio
figure(2)
plot(potention_ra12,total_infected,'LineWidth',2)
set(gca, 'XScale', 'log');
xlabel('ratio of kid:old')
ylabel('total infected at the end')
title('total infected vs ratio of kid:old')

% loss of expected life years and ratio
% checked from Provisional life expectancy
average_ageG1=12.5;
expect_lossG1=64.2;
average_ageG2=74.5;
expect_lossG2=11.5;
average_ageG3=38.7;
expect_lossG3=39.5;
LOL=zeros(length(potention_ra12),1);
for i=1:length(potention_ra12)
    LOL(i)=total_death_G1(i)*expect_lossG1+total_death_G2(i)*expect_lossG2+total_death_G3(i)*expect_lossG3;
end
figure(3)
plot(potention_ra12,LOL,'LineWidth',2)
set(gca, 'XScale', 'log');
xlabel('ratio of kid:old')
ylabel('loss of expected life years')
title('loss of expected life years vs ratio of kid:old')
%end




