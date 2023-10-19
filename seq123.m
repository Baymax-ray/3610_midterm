function [I_G1,I_G2,I_G3,V_I_G1,V_I_G2,V_I_G3,D_G1,D_G2,D_G3,V_D_G1,V_D_G2,V_D_G3,NV_P_G1,NV_P_G2,NV_P_G3]=seq123(VDR,VIR,V,DinM,iniP_G1,iniP_G2,iniP_G3,DR_G1,DR_G2,DR_G3,IR_G1toG1,IR_G1toG2,IR_G1toG3,IR_G2toG1,IR_G2toG2,IR_G2toG3,IR_G3toG1,IR_G3toG2,IR_G3toG3,R_G,FLW)
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
I_G1=[1]; %initial infected population
I_G2=[1];
I_G3=[1];
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
FLW_V=false;

%% time to stimulate
while (I_G1(end)+I_G2(end)+I_G3(end)+V_I_G1(end)+V_I_G2(end)+V_I_G3(end))>=0 && (P_G1+P_G2+P_G3)>0

    %give vaccine first to kids, then old, then mid, assume the vaccine is evenly distributed each day
    %give to FLW first, assume FLW belongs to G3
    D_amount=V/DinM;
    temp1=0;
    temp2=0;
    temp3=0;
    if FLW_V==false
        if D_amount>=FLW
            D_amount=D_amount-FLW;
            temp3=FLW;
            FLW_V=true;
        else
            temp3=D_amount;
            FLW=FLW-D_amount;
            D_amount=0;
        end
    end
    if P_G1>0 && D_amount>0 %give vaccine to g1 first
        if D_amount>P_G1 %we have enough vaccine to give to g1
            D_amount=D_amount-P_G1;
            temp1=P_G1;
        else
            temp1=D_amount;
            D_amount=0;
        end
    end
    if P_G2>0 && D_amount>0
        if D_amount>P_G2 %we have enough vaccine to give to g2
            D_amount=D_amount-P_G2;
            temp2=P_G2;
        else
            temp2=D_amount;
            D_amount=0;
        end
    end
    if P_G3>0 && D_amount>0
        if D_amount>P_G3 %we have enough vaccine to give to g3
            D_amount=D_amount-P_G3;
            temp3=temp3+P_G3;
        else
            temp3=temp3+D_amount;
            D_amount=0;
        end
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

