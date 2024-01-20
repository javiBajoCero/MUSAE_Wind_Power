clear;clc;
vw=7;
rho=1.225;
N=80;
D=76;
poles=2;
Ugrid=960;
f=50;
Pn = 2*10^6;

we=2*pi*f;
ws=we/poles;
wg=ws*(1-ss);
wt=wg/N;
lam= wt*(D/2) / vw;
Cpp= 0.0045 * (100 - (lam-10)^2)
P1=0.5*rho*Cpp*(pi*(D^2)/4).*vw^3

Rs=0.005;
Xs=2*pi*50*4e-4;
Rm=140;
Xm=2*pi*50*15e-3;
Rr=0.009;
Xr=2*pi*50*3e-4;

% Equivalent model
Zs=Rs+j*Xs;Zm=Rm*(j*Xm)/(Rm+j*Xm);Zr=Rr+j*Xr;
Zsm=Zs*Zm/(Zs+Zm);
Zth=Zr+Zsm;
V=Ugrid/sqrt(3);
Uth= Zm*V/(Zm+Zs);
Uu=abs(Uth);

% Initialization
ss=0;kk=0;err=1;P1=0;
while ((kk<100)&&(err>1e-6))
 kk=kk+1;
 we=2*pi*f;ws=we/poles;wg=ws*(1-ss);wt=wg/N;
 lam= wt*(D/2) / vw;
 Cpp= 0.0045 * (100 - (lam-10)^2);
 P1a=min(0.5*rho*Cpp*(pi*(D^2)/4)*vw^3,Pn);
 err=abs(P1a-P1);
 P1=P1a;
 Pmech=-P1/3;
 a1=1;a2=2*real(Zth)*Pmech-Uu^2;a3=abs(Zth)^2*(Pmech^2);
 sol=roots([a1 0 a2 0 a3]);
 Rt=sol.^2/Pmech;
 s=Rr./(Rr+Rt);
 ss=s(2);
end;
P1
ss
Igrid = V/(Zs+Zm*(Zr+Rt(2))/(Zm+(Zr+Rt(2))))
Sgrid= 3*V*conj(Igrid)
losses=real(Sgrid)+P1


ssaux=-1:.00001:1;
Iaux=Uth./(Zth+Rr*(1-ssaux)./ssaux);
Paux=3*abs(Iaux).^2.*Rr.*(1-ssaux)./ssaux;
wgaux=ws*(1-ssaux);
wtaux=wgaux/N;
lamaux= wtaux*(D/2) / vw;
Cppaux= 0.0045 * (76 - (lamaux-8.83).^2);
Pauxmec=min(0.5*rho*Cppaux*(pi*(D^2)/4)*vw^3,Pn);
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





Rs=0.005;
Xs=2*pi*50*4e-4;
Rm=140;
Xm=2*pi*50*15e-3;
Rr=0.009;
Xr=2*pi*50*3e-4;




%% Develop a Matlab program to locate the operating points in the Cplambda curve.
lam1=0:.1:18;
Cpp1= 0.0045 * (100 - (lam1-10).^2);
figure(1);
plot(lam1,Cpp1);hold on;grid on;
plot(lam,Cpp,'ko');
for ii = 1:1:length(vw)
    txt{ii} = ['v_w = ' num2str(vw(ii)) ' m/s'];
    text(lam(ii), Cpp(ii), txt{ii}, 'FontSize', 12);
end;

xlabel('\lambda','FontSize',12);
ylabel('C_p','FontSize',12);

%% Develop a Matlab program to plot the power generated for different wind
%% speeds. Include a comparison between the obtained power and the maximum 
%available power.


figure(2);
vw2=0:.1:15
lam2= min(20,max(wt*(D/2) ./ vw2,0));
Cpp2= 0.0045 * (100 - (lam2-10).^2)
P2=0.5*rho*Cpp2*(pi*(D^2)/4).*vw2.^3;
P2a=min(0.5*rho*Cpp2*(pi*(D^2)/4).*vw2.^3,Pn);
h=subplot(1,1,1);
plot(vw2,P2,':k');hold on;grid on;
plot(vw2,P2a,'k');
plot(vw,P1,'ko');
for ii=1:1:1
 txt{ii}=['v_w = ' num2str(vw(ii)) ' m/s'];
 text(vw(ii),P1(ii),txt{ii},'FontSize',18);
end;
xlabel('v_w','FontSize',18);
ylabel('P','FontSize',18);
set(h,'FontSize',18);


figure(3);
lam2b= 10;
Cpp2b= 0.45;
P3=0.5*rho*Cpp2b*(pi*(D^2)/4).*vw2.^3;
P3a=min(0.5*rho*Cpp2b*(pi*(D^2)/4).*vw2.^3,Pn);
h=subplot(1,1,1);
plot(vw2,P2,':k');hold on;grid on;
plot(vw2,P2a,'k','LineWidth',3);
plot(vw,P1,'ko');
plot(vw2,P3,':r');
plot(vw2,P3a,'r');
for ii=1:1:1
 txt{ii}=['v_w = ' num2str(vw(ii)) ' m/s'];
 text(vw(ii),P1(ii),txt{ii},'FontSize',18);
end;
xlabel('v_w','FontSize',18);
ylabel('P','FontSize',18);
set(h,'FontSize',18);
axis([3 15 0 6e6]);

%%new cp
P = real(Sgrid);  % Potencia activa
Q = imag(Sgrid);  % Potencia reactiva

Eff=P*100/P1;

% Mostrar los resultados
disp(['Slip: ' num2str(ss) ' slip']);
disp(['CP: ' num2str(Cpp) ' cp']);
disp(['Pmech: ' num2str(P1) ' W']);
disp(['Tip speed ratio: ' num2str(lam) ' lambda']);
disp(['Wgen: ' num2str(wg) ' rad/s']);
disp(['Wt: ' num2str(wt) ' rad/s']);
disp(['Potencia Activa (P): ' num2str(-P) ' Watts']);
disp(['Potencia Reactiva (Q): ' num2str(-Q) ' VAR']);
disp(['Generator eff: ' num2str(-Eff) ' %']);