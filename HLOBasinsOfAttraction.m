clear all
%HLO Basins of Attraction

global Ts b a d f m h s g u v B

% Bistable HLO EP/EP Paramaters (0.7585, 0, 0.1676) and (0.0233, 0, 0.1122)
B=0.0025; %(0.0049,0,09); 0.0625 in Scenario 1
a=0.2; %(0.4,2.7); a in Scenario 1
d=0.08; %(0.004,0.6); 0.08 in Scenario 1
h=0.07; %(0.33,1) maybe equal to g; 0.07 in Scenario 1
u=0.5; %(0.07,0.9); 0.5 in Scenario 1
g=0.07; %(0.33,1) maybe equal to h; 0.07 in Scenario 1
v=1/3; %v>0 in scenario 1; 1/3 in Scenario 1
b=2;
f=1;
m=0.2;
s=1/3;
Ts=0.65; %0.47 red 0.48, 0.49 green, 0.5 cyan

% %Bistable HLO Cycle/EP Paramaters (0.146, 0, 0.116) and (0.035, 0, 0.032)
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

Tmax=1000;
L=0.1;
Tol=0.001;

Hare=[L:L:1];
Lynx=[L:L:1];
Owl=[L:L:1];

HSS=m*b/(f-m);
OSS=-Ts*s/(Ts*h*HSS^2/(B+HSS^2)+(1-Ts)*g*HSS/(d+HSS)-u)-v;
LSS=(b+HSS)*(Ts*(1-HSS-HSS*OSS/(B+HSS^2))-(1-Ts)*a*OSS/(d+HSS));

figure(1)
xlabel('Hare')
ylabel('Owl')
clear title
titlestr=['Basins of Attraction'];
title(titlestr)
hold on

for i=1:length(Lynx)
    for j=1:length(Owl)
        
        clear X T
        x0=Hare(i);
        y0=0.01 ;
        z0=Owl(j);
        [T,X]=ode45(@HLOode,0:0.001:Tmax,[x0,y0,z0]);
        
        if abs(HSS-X(end,1))<=Tol && abs(LSS-X(end,2))<=Tol && abs(OSS-X(end,3))<=Tol
            figure(1)
            plot(x0,z0,'b.')
            hold on
        else
            figure(1)
            plot(x0,z0,'g.')
            hold on
        end
    end
end
        