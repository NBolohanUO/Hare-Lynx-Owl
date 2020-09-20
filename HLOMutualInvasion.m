clear GLOBAL
%HLO Mutual Invasion

global b B a d f m h v u g

B=0.0625; %(0.0049,0,09)
a=0.2; %(0.4,2.7)
d=0.08; %(0.004,0.6)
h=0.07; %(0.33,1) maybe equal to g
u=0.5; %(0.07,0.9)
g=0.07; %(0.33,1) maybe equal to h
v=1/3; %v>0 in scenario 1
b=0.3;
f=3.2;
m=1.5;

sres=0.74;

Jump=0.01; %Grid step size
RosMacCond=f-m-b*(f+m); %>0 for stable steady state in xy system, <0 for periodic orbit

Tinit=500;
options = odeset('RelTol',1e-8,'AbsTol',[1e-6 1e-6 1e-6]);

Tmax=500;

PC=0; %Percentage trackers
PCold=0;

[TT,SS]=meshgrid(Jump:Jump:1-Jump,sres*Jump:sres*Jump:sres);
[J,I]=size(TT);
MU=zeros(size(TT)); %Floquet multiplier
OIC=zeros(size(TT));

for i=1:I
    for j=1:J
        Ts=TT(1,i);
        s=SS(j,1);
        
        clear X T
        
        x0=0.5;
        y0=0;
        z0=0.5;
        
        [T,X]=ode45(@HLOode,[0 Tmax],[x0 y0 z0],options);
        xbar=X(end,1);
        zbar=X(end,3);
        LIC(j,i)=f*xbar/(b+xbar)-m; 
        PCold=PC;
        PC=floor(100*((i-1)*J+j)/(I*J));
        if PC==PCold+1
            disp(['Lynx InvCond (xbar=',num2str(xbar),', RMC=',num2str(RosMacCond),'): ',num2str(PC),'% done.'])
        end
    end
end

LIC=LIC;
figure(1)
contour(TT,SS,LIC,[0 0],'m')
% hold on
% contour(TT,SS,LIC,[0.3 0.3],'m--')
hold on
axis([0 1 0 sres])

PC=0; %Percentage trackers
PCold=0;

Tinit=500; %If RosMacCond<0, run xy0 for Tinit to get on orbit, then extract orbit during Tfinal time
Tfinal=100;

if RosMacCond<0
    xbar=m*b/(f-m); %xbar value in xy system
    disp(['xbar=',num2str(xbar),'. RM Condition = ',num2str(RosMacCond)])
    OIC=TT.*(h.*xbar^2./(B+xbar^2)+SS./v-u)+(1-TT).*(g.*xbar./(d+xbar)-u); %Owl Invasion Condition
    figure(1)
    xlabel('T')
    ylabel('s')
    clear title
    titlestr=['Invasion Conditions (fp). B=',num2str(B),', b=',num2str(b),', d=',num2str(d),', g=',num2str(g),', h=',num2str(h),', m=',num2str(m),', u=',num2str(u),', v=',num2str(v)];
    title(titlestr)
    hold on
    contour(TT,SS,OIC,[0 0],'g')
%     hold on
%     contour(TT,SS,OIC,[0.1 0.1],'g--')
    hold on
    axis([0 1 0 sres])
elseif RosMacCond>0
    IndexPeaks=0;
    
    for i=1:I
        for j=1:J
            Ts=TT(1,i);
            s=SS(j,1);
            
            x0=0.5;
            y0=0.5;
            z0=0;
            
            clear X T
            
            [T,X]=ode45(@HLOode,[0 Tinit],[x0,y0,z0],options);
            
            x0 = X(end,1); %Get on periodic orbit and take endpoint
            y0 = X(end,2);
            z0 = X(end,3);

            while length(IndexPeaks)<2 %Run system until we get >1 full orbit. Increase Tfinal until we do
                Tfinal=Tfinal+100;
                clear X T
            
                [T,X]=ode45(@HLOode,[0 Tfinal],[x0,y0,z0],options);
            
                Tinter=linspace(T(1),T(end),100*Tfinal);
                Hinter=interp1(T,X(:,1),Tinter);
                Linter=interp1(T,X(:,2),Tinter);
            
                [HPeaks,IndexPeaks]=findpeaks(Hinter);
            end
            
            PerStart=Tinter(IndexPeaks(1));
            PerEnd=Tinter(IndexPeaks(2));
            Period=PerEnd-PerStart;
            
            CycleT=Tinter(IndexPeaks(1):IndexPeaks(2));
            CycleH=Hinter(IndexPeaks(1):IndexPeaks(2));
            CycleL=Linter(IndexPeaks(1):IndexPeaks(2));
            
            A=Ts.*(h.*CycleH.^2./(B+CycleH.^2)+s./v-u)+(1-Ts)*(g*CycleH./(d+CycleH)-u);
            MU(j,i)=(1/Period)*trapz(A);
            
            PCold=PC;
            PC=floor(100*((i-1)*J+j)/(I*J));
            if PC==PCold+1
                disp(['Owl InvCond (Periodic Orbit, xbar=',num2str(xbar),', RMC=',num2str(RosMacCond),'): ',num2str(PC),'% done.'])
            end
        end
    end
    
    figure(1)
    xlabel('T')
    ylabel('s')
    clear title
    titlestr=['Invasion Conditions -(po). B=',num2str(B),', b=',num2str(b),', d=',num2str(d),', g=',num2str(g),', h=',num2str(h),', m=',num2str(m),', u=',num2str(u),', v=',num2str(v)];
    title(titlestr)
    hold on
    contour(TT,SS,MU,[0 0],'g')
%     hold on
%     contour(TT,SS,MU,[10 10],'g--')
    hold on
    axis([0 1 0 sres])
else
    disp('RosMacCond=0, zero eigenvalue in Hare-Lynx system.')
end        