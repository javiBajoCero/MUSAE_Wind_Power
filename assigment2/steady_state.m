
%The rationale of Problem 1 can be followed with the only difference of not assuming s=0.
%now the ws is not equal to wg   ( we|ws|wg|wt ) we take intoa account the
%slip


clc;
clear;
close all;
% Parameters
rho=1.225;
N=80; %transmission ratio?
D=76; %blade sweep diameter
poles=2;


Ugrid=960;
f=50; omega=2*pi*f;
vw=14;

Rs=0.005;
Xs=omega*4e-4;
Rr=0.009;
Xr=omega*3e-4;
Rm=140;
Xm=omega*15e-3;

% Equivalent model
Zs=Rs+1i*Xs;Zm=Rm*(1i*Xm)/(Rm+1i*Xm);Zr=Rr+1i*Xr;
Zsm=Zs*Zm/(Zs+Zm);
Zth=Zr+Zsm;
V=Ugrid/sqrt(3);
Uth= Zm*V/(Zm+Zs);
Uu=abs(Uth);



% Initialization, we are solving numerically because solving the ecuations
% are tough
ss=0; %0 power
kk=0;
err=1;
P1=0;
while ((kk<100)&&(err>1e-6))
 kk=kk+1;

 we=2*pi*f;   %chain of speeds, this is how they relate
 ws=we/poles;
 wg=ws*(1-ss);
 wt=wg/N;

 lam= wt*(D/2) / vw;    %ASSUMING WT SO WE CAN COMPUTE A LAMBDA

 Cpp= 0.0045 * (100 - (lam-10)^2);  %with lambda we calculate cpp power coeficient
 P1a=0.5*rho*Cpp*(pi*(D^2)/4)*vw^3; %power

 err=abs(P1a-P1);                   %this establishes the end of the while, compares past power with current power
 P1=P1a;
 %now i calculate the new slip
 Pmech=-P1/3; %power mechanical
 a1=1;
 a2=2*real(Zth)*Pmech-Uu^2;
 a3=abs(Zth)^2*(Pmech^2);
 sol=roots([a1 0 a2 0 a3]); %U
 Rt=sol.^2/Pmech;
 s=Rr./(Rr+Rt); %computing new slip
 ss=s(2);
end

P1
ss
Igrid = V/(Zs+Zm*(Zr+Rt(2))/(Zm+(Zr+Rt(2))))
Sgrid= 3*V*conj(Igrid)
losses=real(Sgrid)+P1
Cpp

ssaux=-1:.00001:1;
Iaux=Uth./(Zth+Rr*(1-ssaux)./ssaux);
Paux=3*abs(Iaux).^2.*Rr.*(1-ssaux)./ssaux;
wgaux=ws*(1-ssaux);
wtaux=wgaux/N;
lamaux= wtaux*(D/2) / vw;
Cppaux= 0.0045 * (100 - (lamaux-10).^2);
Pauxmec=0.5*rho*Cppaux*(pi*(D^2)/4)*vw^3;

figure(1);

h=subplot(1,1,1);
plot(wgaux,Paux,'k');hold on;grid on;
plot(wgaux,-Pauxmec,'r');
plot(wg,-P1,'ok');
xlabel('\omega [rad/s]','FontSize',18);
ylabel('P [W]','FontSize',18);
legend('Electrical power','Mechanical power','Operating point');
set(h,'FontSize',18);
axis([0 200 -3e6 3e6]);