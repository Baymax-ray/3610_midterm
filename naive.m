function [I_G1,I_G2,I_G3,V_I_G1,V_I_G2,V_I_G3,D_G1,D_G2,D_G3,V_D_G1,V_D_G2,V_D_G3,NV_P_G1,NV_P_G2,NV_P_G3]=naive(VDR,VIR,V,DinM,iniP_G1,iniP_G2,iniP_G3,DR_G1,DR_G2,DR_G3,IR_G1toG1,IR_G1toG2,IR_G1toG3,IR_G2toG1,IR_G2toG2,IR_G2toG3,IR_G3toG1,IR_G3toG2,IR_G3toG3,R_G)
% people always acccept vaccine
% no priority

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
% total vaccinated population
NV_P_G1=[0];
NV_P_G2=[0];
NV_P_G3=[0];

%% time to stimulate
while (I_G1(end)+I_G2(end)+I_G3(end)+V_I_G1(end)+V_I_G2(end)+V_I_G3(end))>=0 && (P_G1+P_G2+P_G3)>0

    %give vaccine proportion to each group, assume the vaccine is evenly distributed each day
    temp=P_G1+P_G2+P_G3;
    temp1=(V/DinM)*P_G1/temp;
    temp2=(V/DinM)*P_G2/temp;
    temp3=(V/DinM)*P_G3/temp;
    if temp1>P_G1
        temp1=P_G1;
    end
    if temp2>P_G2
        temp2=P_G2;
    end
    if temp3>P_G3
        temp3=P_G3;
    end

    V_P_G1=V_P_G1+temp1;
    V_P_G2=V_P_G2+temp2;
    V_P_G3=V_P_G3+temp3;
    NV_P_G1=[NV_P_G1 temp1];
    NV_P_G2=[NV_P_G2 temp2];
    NV_P_G3=[NV_P_G3 temp3];

    P_G1=P_G1-temp1;
    P_G2=P_G2-temp2;
    P_G3=P_G3-temp3;

    %get infected
    startIndex = max(1, numel(I_G1) - R_G+1);
    
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
end
