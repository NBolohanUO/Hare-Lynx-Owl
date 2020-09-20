clear all
%HLO Owl Simp Steady State Equations and Eigenvalues

global a b k d f m h g u

%Steady State
k=0.5;
d=1;
a=0.5;
b=0.8;
f=2;
m=1;
h=0.5;
g=1;
u=1.8;

%Periodic Orbit
k=0.9;
d=1;
a=1;
b=0.5;
f=2;
m=0.3;
h=0.3;
g=1;
u=0.2;

Jump=0.01;
L=0;

TT=[Jump:Jump:1];
SS=[5*Jump:5*Jump:5];

figure(1)
xlabel('T')
ylabel('s')
% clear title
% titlestr=['Steady State Equations (Owl Simp). k=',num2str(k),', d=',num2str(d),' a=',num2str(a),', b=',num2str(b),', f=',num2str(f),', m=',num2str(m),', h=',num2str(h),', g=',num2str(g),', u=',num2str(u)];
% title(titlestr)
hold on

figure(2)
xlabel('T')
ylabel('s')
% clear title
% titlestr=['Eigenvalues (Owl Simp). k=',num2str(k),', d=',num2str(d),', a=',num2str(a),', b=',num2str(b),', f=',num2str(f),', m=',num2str(m),', h=',num2str(h),', g=',num2str(g),', u=',num2str(u)];
% title(titlestr)
hold on

PC=0;
   
for i=1:length(TT)
    for j=1:length(SS)
        T=TT(i);
        s=SS(j);
        
        x=m*b/(f-m*k);
        z=(T*h*x^2+(1-T)*g*x+T*s-(1-T)*u)/(T*s);
        y=(b+k*x)*(T-T*x-T*a*x*z+(1-T)*(-d*z));
        
        if x>L && y>L && z>L %HLO, black
            figure(1)
            plot(T,s,'k.')
            hold on
            
        elseif x>L && y>L && z<=L %HL, magenta
            figure(1)
            plot(T,s,'m.')
            hold on
            
        elseif x>L && y<=L && z>L %HO, green
            figure(1)
            plot(T,s,'g.')
            hold on
            
        elseif x<=L && y>L && z>L %LO, cyan
            figure(1)
            plot(T,s,'c.')
            hold on
            
        elseif x>L && y<=L && z<=L %H, blue
            figure(1)
            plot(T,s,'b.')
            hold on
            
        elseif x<=L && y>L && z<=L %L, red
            figure(1)
            plot(T,s,'r.')
            hold on
            
        elseif x<=L && y<=L && z>L %O, yellow
            figure(1)
            plot(T,s,'y.')
            hold on
            
        else x<=L && y<=L && z<= L %None, white
            figure(1)
            plot(T,s,'w.')
            hold on
        end
        
        J=[T*(1-2*x-2*a*x*z)+(1-T)*(-d*z)-b*y/(b+k*x)^2 -x/(b+k*x) T*(-a*x^2)+(1-T)*(-d*x);
            f*b*y/(b+k*x)^2 f*x/(b+k*x)-m 0;
            T*(2*h*x*z)+(1-T)*(g*z) 0 T*(h*x^2+s-2*s*z)+(1-T)*(g*x-u)];
       
       E=eig(J);
       e1=E(1);
       e2=E(2);
       e3=E(3);
       
       if real(e1)<0 && real(e2)<0 && real(e3)<0 %All negative, black
           if any(imag(E)~=0)
               figure(2)
               plot(T,s,'k+')
               hold on
           else
               figure(2)
               plot(T,s,'k.')
               hold on
           end
           
       elseif (real(e1)<0 && real(e2)<0 && real(e3)>=0)||(real(e1)<0 && real(e2)>=0 && real(e3)<0)||(real(e1)>=0 && real(e2)<0 && real(e3)<0) %One positive, red
           if any(imag(E)~=0)
               figure(2)
               plot(T,s,'r+')
               hold on
           else
               figure(2)
               plot(T,s,'r.')
               hold on
           end
       
       elseif (real(e1)<0 && real(e2)>=0 && real(e3)>=0)||(real(e1)>=0 && real(e2)<0 && real(e3)>=0)||(real(e1)>=0 && real(e2)>=0 && real(e3)<0) %Two positive, cyan
           if any(imag(E)~=0)
               figure(2)
               plot(T,s,'c+')
               hold on
           else
               figure(2)
               plot(T,s,'c.')
               hold on
           end
       
       else %All positive, green
           if any(imag(E)~=0)
               figure(2)
               plot(T,s,'g+')
               hold on
           else
               figure(2)
               plot(T,s,'g.')
               hold on
           end
       end
       
       PCold=PC;
        PC=floor(100*((i-1)*length(SS)+j)/(length(TT)*length(SS)));
        if PC==PCold+1
            disp(['Steady State Eqns and Evals: ',num2str(PC),'% done.'])
        end
        
    end
end

% zstar= @(T,s) (h/s)*((m/f)^2)+1+((1-T)/(T*s))*(g*(m/f)-u);
% ystar= @(T,s) (T/c)*(1-(m/f)-a*(m/f)*(h/s)*((m/f)^2)+1+((1-T)/(T*s))*(g*(m/f)-u))+((1-T)/c)*(-b*(h/s)*((m/f)^2)+1+((1-T)/(T*s))*(g*(m/f)-u));



% Contour Plots for y* and z*

% figure(2)
% fimplicit(zstar,[0 1 0 5], 'g')
% hold on
% fimplicit(ystar,[0 1 0 5], 'r')
% hold on
% xlabel('T')
% ylabel('s')
% clear title
% titlestr=['a=',num2str(a),', c=',num2str(c),', b=',num2str(b),', f=',num2str(f),', m=',num2str(m),', h=',num2str(h),', g=',num2str(g),', u=',num2str(u)];
% title(titlestr)

% Jump=0.01; %Grid step size
% xbar=m*b/(f-m*k); %xbar value in xy system
% RosMacCond=f*(k-b)-m*k*(k+b); %<0 for stable steady state in xy system, >0 for periodic orbit
% 
% Tinit=3000; %If RosMacCond>0, run xy0 for Tinit to get on orbit, then extract orbit during Tfinal time
% Tfinal=300;
% 
% PC=0; %Percentage trackers
% PCold=0;
% 
% [TT,SS]=meshgrid(Jump:Jump:1-Jump,5*Jump:5*Jump:10);
% [J,I]=size(TT);
% MU=zeros(size(TT)); %Floquet multiplier
% 
% % Lynx Invasion Stuff
% 
% C0=@(s,T) T.*(1-T).*s.*d-(1-T).^2.*u.*d-T.^2.*s;
% C1=@(s,T) T.^2.*s.*a+(1-T).^2.*g.*d-T.*(1-T).*u.*a+T.^2.*s;
% C2=@(s,T) T.*(1-T).*(h.*d+g.*a);
% C3=@(s,T) T.^2.*h.*a;
% 
% R=NaN(size(TT)); %xbar in x0z system
% 
% HOCond1=zeros(size(TT));
% HOCond2=zeros(size(TT));
% 
% for i=1:I
%     for j=1:J
%         Ts=TT(1,i);
%         s=SS(j,1);
%         
%         HOCond1(j,i)=Ts/((1-Ts)*d)-(Ts*s-(1-Ts)*u)/(Ts*s);
%         HOCond2(j,i)=(Ts*(h+s)+(1-Ts)*(g-u))/(Ts*s);
%         
%         if (HOCond1(j,i)>0)&&(HOCond2(j,i)>0)
%             r=roots([C3(s,Ts) C2(s,Ts) C1(s,Ts) C0(s,Ts)]);
%             if any(imag(r)~=0)
%                 R(j,i)=r(find(imag(r)==0));
%             else
%                 LessOne=real(r)<=1;
%                 BigZero=real(r)>=0;
%                 R(j,i)=r(find(LessOne.*BigZero==1));
%             end
%         end
%     end
% end
% 
% LynxInvCond=R-xbar*ones(size(R));
% 
% figure(1)
% contour(TT,SS,LynxInvCond,[0 0],'r')
% hold on
% % contour(TT,SS,LynxInvCond,[0.01 0.01],'m--')
% % hold on
% % contour(TT,SS,HOCond1,[0 0],'b')
% % hold on
% % % contour(TT,SS,HOCond1,[1 1],'b--')
% % % hold on
% % contour(TT,SS,HOCond2,[0 0],'r')
% % hold on
% % % contour(TT,SS,HOCond2,[0.1 0.1],'r--')
% % % hold on
% axis([0 1 0 5])
% 
% if RosMacCond<0
%     disp(['xbar=',num2str(xbar),'. RM Condition = ',num2str(RosMacCond)])
%     OIC=TT.*(h.*(m.*b./(f-m.*k)).^2+SS)+(1-TT).*(g.*m.*b./(f-m.*k)-u); %Owl Invasion Condition
%     figure(1)
%     xlabel('T')
%     ylabel('s')
% %     clear title
% %     titlestr=['Invasion Conditions (Owl Simp, Steady State). k=',num2str(k),', d=',num2str(d),' a=',num2str(a),', b=',num2str(b),', f=',num2str(f),', m=',num2str(m),', h=',num2str(h),', g=',num2str(g),', u=',num2str(u)];
% %     title(titlestr)
%     hold on
%     contour(TT,SS,OIC,[0 0],'b')
%     hold on
% %     contour(TT,SS,OIC,[10 10],'g--')
% %     hold on
% elseif RosMacCond>0
%     IndexPeaks=0;
%     
%     options = odeset('RelTol',1e-8,'AbsTol',[1e-6 1e-6 1e-6]);
%     
%     for i=1:I
%         for j=1:J
%             Ts=TT(1,i);
%             s=SS(j,1);
%             
%             x0=0.5;
%             y0=0.5;
%             z0=0;
%             
%             clear X T
%             
%             [T,X]=ode45(@(t,x) [Ts*(x(1)*(1-x(1))-x(1)*x(2)/(b+k*x(1))-a*x(1)^2*x(3))+(1-Ts)*(-x(1)*x(2)/(b+k*x(1))-d*x(1)*x(3)); f*x(1)*x(2)/(b+k*x(1))-m*x(2); Ts*(h*x(1)^2*x(3)+s*x(3)*(1-x(3)))+(1-Ts)*(g*x(1)*x(3)-u*x(3))],[0 Tinit],[x0,y0,z0],options);
%             
%             x0 = X(end,1); %Get on periodic orbit and take endpoint
%             y0 = X(end,2);
%             z0 = X(end,3);
%             Tfinal=300;
%             
%             while length(IndexPeaks)<2 %Run system until we get >1 full orbit. Increase Tfinal until we do
%                 Tfinal=Tfinal+100;
%                 clear X T
%             
%                 [T,X]=ode45(@(t,x) [Ts*(x(1)*(1-x(1))-x(1)*x(2)/(b+k*x(1))-a*x(1)^2*x(3))+(1-Ts)*(-x(1)*x(2)/(b+k*x(1))-d*x(1)*x(3)); f*x(1)*x(2)/(b+k*x(1))-m*x(2); Ts*(h*x(1)^2*x(3)+s*x(3)*(1-x(3)))+(1-Ts)*(g*x(1)*x(3)-u*x(3))],[0 Tfinal],[x0,y0,z0],options);
%             
%                 Tinter=linspace(T(1),T(end),100*Tfinal);
%                 Hinter=interp1(T,X(:,1),Tinter);
%                 Linter=interp1(T,X(:,2),Tinter);
%             
%                 [HPeaks,IndexPeaks]=findpeaks(Hinter);
%             end
%             
%             PerStart=Tinter(IndexPeaks(1));
%             PerEnd=Tinter(IndexPeaks(2));
%             Period=PerEnd-PerStart;
%             
%             CycleT=Tinter(IndexPeaks(1):IndexPeaks(2));
%             CycleH=Hinter(IndexPeaks(1):IndexPeaks(2));
%             CycleL=Linter(IndexPeaks(1):IndexPeaks(2));
%             
%             A=Ts*(h*CycleH.^2+s)+(1-Ts)*(g*CycleH-u);
%             MU(j,i)=(1/Period)*trapz(A);
%             
%             PCold=PC;
%             PC=floor(100*((i-1)*J+j)/(I*J));
%             if PC==PCold+1
%                 disp(['Owl InvCond (Periodic Orbit, xbar=',num2str(xbar),', RMC=',num2str(RosMacCond),'): ',num2str(PC),'% done.'])
%             end
%         end
%     end
%     
%     figure(1)
%     xlabel('T')
%     ylabel('s')
% %     clear title
% %     titlestr=['Invasion Conditions (Owl Simp, Orbit). k=',num2str(k),', d=',num2str(d),' a=',num2str(a),', b=',num2str(b),', f=',num2str(f),', m=',num2str(m),', h=',num2str(h),', g=',num2str(g),', u=',num2str(u)];
% %     title(titlestr)
%     hold on
%     contour(TT,SS,MU,[0 0],'b')
%     hold on
% %     contour(TT,SS,MU,[3 3],'g--')
% %     hold on
% else
%     disp('RosMacCond=0, zero eigenvalue in Hare-Lynx system.')
% end





