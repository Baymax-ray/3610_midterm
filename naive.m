%% initial parameters
P_tot=107000;
V=3500;
DinM=30;% 30 days in a month
iniP_G1=P_tot*0.137; %<18
iniP_G2=P_tot*0.164; %>65
iniP_G3=P_tot-iniP_G1-iniP_G2; %18-65
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
%death rate
DR_G1=0.001;
DR_G2=0.18;
DR_G3=0.02;
%infection rate
IR_G1toG1=0.18; 
IR_G1toG2=0.28; 
IR_G1toG3=0.18; 
IR_G2toG1=0.17; 
IR_G2toG2=0.25; 
IR_G2toG3=0.17; 
IR_G3toG1=0.08; %the ratio of parents to children : children to parents is about 3:7
IR_G3toG2=0.28; 
IR_G3toG3=0.18; 

VDR=0.1; %vaccine effect on death rate
VIR=0.3; %vaccine effect on infection rate

%recover days
R_G=10; %assume everyone is no longer contagious after 10 days

%% time to stimulate
while (I_G1(end)+I_G2(end)+I_G3(end)+V_I_G1(end)+V_I_G2(end)+V_I_G3(end))>=0 && (P_G1+P_G2+P_G3)>0

    %give vaccine proportion to each group, assume the vaccine is evenly distributed each day
    temp=P_G1+P_G2+P_G3;

    V_P_G1=V_P_G1+(V/DinM)*P_G1/temp;
    V_P_G2=V_P_G2+(V/DinM)*P_G2/temp;
    V_P_G3=V_P_G3+(V/DinM)*P_G3/temp;

    P_G1=P_G1-(V/DinM)*P_G1/temp;
    if P_G1<0
        P_G1=0;
    end
    P_G2=P_G2-(V/DinM)*P_G2/temp;
    if P_G2<0
        P_G2=0;
    end
    P_G3=P_G3-(V/DinM)*P_G3/temp;
    if P_G3<0
        P_G3=0;
    end

    %get infected
    startIndex = max(1, numel(I_G1) - R_G+1);
    
    %sumLast10 = sum(arr(startIndex:end));
    new_I_G1=(P_G1/iniP_G1)*(IR_G1toG1*sum(I_G1(startIndex:end))+IR_G2toG1*sum(I_G2(startIndex:end))+IR_G3toG1*sum(I_G3(startIndex:end))+VIR*(sum(V_I_G1(startIndex:end))*IR_G1toG1+sum(V_I_G2(startIndex:end))*IR_G2toG1+sum(V_I_G3(startIndex:end))*IR_G3toG1));
    new_I_G2=(P_G2/iniP_G2)*(IR_G1toG2*sum(I_G1(startIndex:end))+IR_G2toG2*sum(I_G2(startIndex:end))+IR_G3toG2*sum(I_G3(startIndex:end))+VIR*(sum(V_I_G1(startIndex:end))*IR_G1toG2+sum(V_I_G2(startIndex:end))*IR_G2toG2+sum(V_I_G3(startIndex:end))*IR_G3toG2));
    new_I_G3=(P_G3/iniP_G3)*(IR_G1toG3*sum(I_G1(startIndex:end))+IR_G2toG3*sum(I_G2(startIndex:end))+IR_G3toG3*sum(I_G3(startIndex:end))+VIR*(sum(V_I_G1(startIndex:end))*IR_G1toG3+sum(V_I_G2(startIndex:end))*IR_G2toG3+sum(V_I_G3(startIndex:end))*IR_G3toG3));

    new_V_I_G1=VIR*(V_P_G1/iniP_G1)*(IR_G1toG1*sum(I_G1(startIndex:end))+IR_G2toG1*sum(I_G2(startIndex:end))+IR_G3toG1*sum(I_G3(startIndex:end))+VIR*(sum(V_I_G1(startIndex:end))*IR_G1toG1+sum(V_I_G2(startIndex:end))*IR_G2toG1+sum(V_I_G3(startIndex:end))*IR_G3toG1));
    new_V_I_G2=VIR*(V_P_G2/iniP_G2)*(IR_G1toG2*sum(I_G1(startIndex:end))+IR_G2toG2*sum(I_G2(startIndex:end))+IR_G3toG2*sum(I_G3(startIndex:end))+VIR*(sum(V_I_G1(startIndex:end))*IR_G1toG2+sum(V_I_G2(startIndex:end))*IR_G2toG2+sum(V_I_G3(startIndex:end))*IR_G3toG2));
    new_V_I_G3=VIR*(V_P_G3/iniP_G3)*(IR_G1toG3*sum(I_G1(startIndex:end))+IR_G2toG3*sum(I_G2(startIndex:end))+IR_G3toG3*sum(I_G3(startIndex:end))+VIR*(sum(V_I_G1(startIndex:end))*IR_G1toG3+sum(V_I_G2(startIndex:end))*IR_G2toG3+sum(V_I_G3(startIndex:end))*IR_G3toG3));

    if new_I_G1>P_G1
        new_I_G1=P_G1;
    end
    P_G1=P_G1-new_I_G1;
    if new_I_G2>P_G2
        new_I_G2=P_G2;
    end
    P_G2=P_G2-new_I_G2;
    if new_I_G3>P_G3
        new_I_G3=P_G3;
    end
    P_G3=P_G3-new_I_G3;

    if new_V_I_G1>V_P_G1
        new_V_I_G1=V_P_G1;
    end
    V_P_G1=V_P_G1-new_V_I_G1;
    if new_V_I_G2>V_P_G2
        new_V_I_G2=V_P_G2;
    end
    V_P_G2=V_P_G2-new_V_I_G2;
    if new_V_I_G3>V_P_G3
        new_V_I_G3=V_P_G3;
    end
    V_P_G3=V_P_G3-new_V_I_G3;

    I_G1=[I_G1 new_I_G1];
    I_G2=[I_G2 new_I_G2];
    I_G3=[I_G3 new_I_G3];

    V_I_G1=[V_I_G1 new_V_I_G1];
    V_I_G2=[V_I_G2 new_V_I_G2];
    V_I_G3=[V_I_G3 new_V_I_G3];

    %get dead, assume all death happend on the first day of infection
    D_G1=[D_G1 DR_G1*I_G1(end)];
    D_G2=[D_G2 DR_G2*I_G2(end)];
    D_G3=[D_G3 DR_G3*I_G3(end)];

    V_D_G1=[V_D_G1 VDR*DR_G1*V_I_G1(end)];
    V_D_G2=[V_D_G2 VDR*DR_G2*V_I_G2(end)];
    V_D_G3=[V_D_G3 VDR*DR_G3*V_I_G3(end)];
end

%% plot result
% plot infected population per day
figure(1)
plot(I_G1+I_G2+I_G3+V_I_G1+V_I_G2+V_I_G3,"Color",[83 81 84]./255)
hold on
plot(I_G1+V_I_G1,"Color",[57 106 177]./255)
plot(I_G2+V_I_G2,"Color",[204 37 41]./255)
plot(I_G3+V_I_G3,"Color",[62 150 81]./255)
xlabel('days')
ylabel('population')
title('new infected population at each day')
legend('total infected population','infected population in group 1','infected population in group 2','infected population in group 3')
hold off

% plot total dead population per day
figure(2)
plot(D_G1+D_G2+D_G3+V_D_G1+V_D_G2+V_D_G3,"Color",[83 81 84]./255)
hold on
plot(D_G1+V_D_G1,"Color",[57 106 177]./255)
plot(D_G2+V_D_G2,"Color",[204 37 41]./255)
plot(D_G3+V_D_G3,"Color",[62 150 81]./255)
xlabel('days')
ylabel('population')
title('new dead population at each day')
legend('total dead population','dead population in group 1','dead population in group 2','dead population in group 3')
hold off

% plot accumulated infected population
figure(3)
plot(cumsum(I_G1+I_G2+I_G3+V_I_G1+V_I_G2+V_I_G3),"Color",[83 81 84]./255)
hold on
plot(cumsum(I_G1+V_I_G1),"Color",[57 106 177]./255)
plot(cumsum(I_G2+V_I_G2),"Color",[204 37 41]./255)
plot(cumsum(I_G3+V_I_G3),"Color",[62 150 81]./255)
xlabel('days')
ylabel('population')
title('accumulated infected population')
legend('total infected population','infected population in group 1','infected population in group 2','infected population in group 3')
hold off

% plot accumulated dead population
figure(4)
plot(cumsum(D_G1+D_G2+D_G3+V_D_G1+V_D_G2+V_D_G3),"Color",[83 81 84]./255)
hold on
plot(cumsum(D_G1+V_D_G1),"Color",[57 106 177]./255)
plot(cumsum(D_G2+V_D_G2),"Color",[204 37 41]./255)
plot(cumsum(D_G3+V_D_G3),"Color",[62 150 81]./255)
xlabel('days')
ylabel('population')
title('accumulated dead population')
legend('total dead population','dead population in group 1','dead population in group 2','dead population in group 3')
hold off
