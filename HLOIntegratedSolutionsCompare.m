clear all
%HLO Integrated Solutions

global Ts b a d f m h s g u v B

%Bistable HLO EP/EP Paramaters
B=0.025; %(0.0049,0,09); 0.0625 in Scenario 1
a=0.2; %(0.4,2.7); a in Scenario 1
d=0.08; %(0.004,0.6); 0.08 in Scenario 1
h=0.07; %(0.33,1) maybe equal to g; 0.07 in Scenario 1
u=0.5; %(0.07,0.9); 0.5 in Scenario 1
g=0.07; %(0.33,1) maybe equal to h; 0.07 in Scenario 1
v=1/3; %v>0 in scenario 1; 1/3 in Scenario 1
b=2;
f=1;
m=0.1753;
s=1/3;
Ts=0.65; %0.47 red 0.48, 0.49 green, 0.5 cyan

% %Bistable HLO Cycle/EP Paramaters
% B=0.0025; %(0.0049,0,09); 0.0625 in Scenario 1
% a=2; %(0.4,2.7); a in Scenario 1
% d=0.1; %(0.004,0.6); 0.08 in Scenario 1
% h=0.05; %(0.33,1) maybe equal to g; 0.07 in Scenario 1
% u=0.4; %(0.07,0.9); 0.5 in Scenario 1
% g=0.1; %(0.33,1) maybe equal to h; 0.07 in Scenario 1
% v=1/3; %v>0 in scenario 1; 1/3 in Scenario 1
% b=2;
% f=1;
% m=0.05;
% s=1/3;
% Ts=0.414; %0.47 red 0.48, 0.49 green, 0.5 cyan

x0=0.4023;
y0=0.2437;
z0=0.153;

Tmax=200;

% figure(1)
% subplot(3,1,1)
% xlabel('Time')
% ylabel('Hare Density')
% subplot(3,1,2)
% xlabel('Time')
% ylabel('Lynx Density')
% subplot(3,1,3)
% xlabel('Time')
% ylabel('Owl Density')
% hold on

figure(1)
xlabel('Time')
ylabel('Hare Density')
hold on

figure(2)
xlabel('Time')
ylabel('Lynx Density')
hold on

figure(3)
xlabel('Time')
ylabel('Owl Density')
hold on

[Tavg,Xavg]=ode15s(@HLOode,0:0.001:Tmax,[x0,y0,z0]);

for i=1:Tmax
    Tsum=Ts;
    Twin=1-Ts;
    
    [T,X]=ode15s(@HLOodeSummer,0:0.001:Tsum,[x0,y0,z0]);
    T=i-1+T;
%     figure(1)
%     subplot(3,1,1)
%     plot(T,X(:,1),'b')
%     hold on
%     subplot(3,1,2)
%     plot(T,X(:,2),'r')
%     hold on
%     subplot(3,1,3)
%     plot(T,X(:,3),'g')
%     hold on

%     figure(1)
%     plot(T,X(:,1),'b')
%     hold on
% 
%     figure(2)
%     plot(T,X(:,2),'r')
%     hold on
% 
%     figure(3)
%     plot(T,X(:,3),'g')
%     hold on

    figure(1)
    plot(i-1+Tsum,mean(X(:,1)),'b*')
    hold on
    
    figure(2)
    plot(i-1+Tsum,mean(X(:,2)),'r*')
    hold on
    
    figure(3)
    plot(i-1+Tsum,mean(X(:,3)),'g*')
    hold on
    
    x0=X(end,1);
    y0=X(end,2);
    z0=X(end,3);
    
    [T,X]=ode15s(@HLOodeWinter,0:0.001:Twin,[x0,y0,z0]);
    T=i-1+Twin+T;
%     figure(1)
%     subplot(3,1,1)
%     plot(T,X(:,1),'b--')
%     hold on
%     subplot(3,1,2)
%     plot(T,X(:,2),'r--')
%     hold on
%     subplot(3,1,3)
%     plot(T,X(:,3),'g--')
%     hold on

%     figure(1)
%     plot(T,X(:,1),'b--')
%     hold on
% 
%     figure(2)
%     plot(T,X(:,2),'r--')
%     hold on
% 
%     figure(3)
%     plot(T,X(:,3),'g--')
%     hold on
    
    figure(1)
    plot(i,mean(X(:,1)),'b*')
    hold on
    
    figure(2)
    plot(i,mean(X(:,2)),'r*')
    hold on
    
    figure(3)
    plot(i,mean(X(:,3)),'g*')
    hold on
    
    x0=X(end,1);
    y0=X(end,2);
    z0=X(end,3);  
end

% figure(1)
% plot(Tavg,Xavg(:,1),'k','LineWidth',2)
% 
% figure(2)
% plot(Tavg,Xavg(:,2),'k','LineWidth',2)
% 
% figure(3)
% plot(Tavg,Xavg(:,3),'k','LineWidth',2)

figure(1)
plot(Tavg,Xavg(:,1),'k','LineWidth',0.5)

figure(2)
plot(Tavg,Xavg(:,2),'k','LineWidth',0.5)

figure(3)
plot(Tavg,Xavg(:,3),'k','LineWidth',0.5)

% figure(1)
% axis([0 Tmax 0 0.7])
% 
% figure(2)
% axis([0 Tmax 0 0.7])
% 
% figure(3)
% axis([0 Tmax 0 0.7])