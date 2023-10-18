%% naive
VDR=0.1;
VIR=0.3;
[I_G1,I_G2,I_G3,V_I_G1,V_I_G2,V_I_G3,D_G1,D_G2,D_G3,V_D_G1,V_D_G2,V_D_G3]=naive(VDR,VIR);




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

% plot total dead population per day
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