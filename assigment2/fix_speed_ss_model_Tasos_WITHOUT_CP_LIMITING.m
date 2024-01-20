clear;
close all;
clc;

% Cp values [ Table ]
c1 = 0.44;
c2 =125;
c3 =0; c4 =0; c5 =0;
c6 = 6.94 ;
c7 = 16.5 ;
c8 =0;
c9 = -0.002;


wind_speed_cut_in   =3;     %[m/s]
wind_speed_cut_off  =20;    %[m/s]
rotor_diameter      =76;    %[m]
rotor_radius        =rotor_diameter/2;%[m]
Area_swept_by_the_blades=pi*rotor_radius^2;    % Area swept by th e b l a d e s
rho=1.225;                  %air density
rotor_nominal_speed =16;    %[1/min]
moment_of_inertia   =9e6;   %[kgm2]
damping_factor      =7.5e6; %[Nm/rad]
transmission_ratio  =80;    %[]
angle_pitch=0;              % blade P i t ch a n g l e



%Electrical parameters of the wind turbine (generator)  
power               =2e6;   %[W]
nominal_voltage     =960;   %[Vph-ph]
grid_freq           =50;    %[Hz]
nominal_current     =1300;  %[A]
connection          ="delta";
pole_pairs          =2;     %[]
rated_speed_50hz    =1500;  %[1/min]
moment_of_inertia   =90;    %[kgm2]
w_ini               =2*pi*grid_freq/(pole_pairs* transmission_ratio) ; % I n i t i a l sp e ed o f th e mach ine
we=2*pi*grid_freq;


stator_resistance   =0.005; %[ohms]
Rs=stator_resistance;
stator_leakage_inductance   =4e-4;%[H]
Xs=we*stator_leakage_inductance;

rotor_resistance    =0.009;%[ohms]
Rr=rotor_resistance;
rotor_leakage_inductance=3e-4;%[H]
Xr=we*stator_leakage_inductance;

magnetizing_inductance=15e-3;%[H]
Xm=we*magnetizing_inductance;
iron_branch_equivalent_resistance=140;%[ohms]
Rm=iron_branch_equivalent_resistance;




vw1=[7 11 14];      %wind speeds
lam= w_ini*(rotor_diameter/2) ./ vw1 %lamda calc
k=(lam+c8*angle_pitch).^(-1)-c9/(1+angle_pitch^3);                             %Auxvariable to calculate Cp
cp=max(0,c1*(c2*k-c3*angle_pitch-c4*angle_pitch^c5-c6).*exp(-c7*k));          %calculation of the Cp for←- different tsrs
Pm=0.5*rho*cp*(pi*(rotor_diameter^2)/4).*vw1.^3      %P mech

Pma=min(0.5*rho*cp*(pi*(rotor_diameter^2)/4).*vw1.^3,power);
cp_new=2.*Pma ./(rho*(pi*(rotor_diameter^2)/4).*vw1.^3);

%% Obtention of the turbine [tipspeedratio(TSR)−Cp] curve
% CP vs tsr curve creation
tsr=0:.1:20;                                                                    %tsr vector from 0 a 17 withincrements of 0.1

for angle_pitch=0:1:29
    pitch=ones(1,29)*angle_pitch;
    weird_constant=1./(tsr+c8.*angle_pitch)-(c9./(1+angle_pitch^3));
    CP=c1.*(c2.*weird_constant -c3.*angle_pitch -c4.*angle_pitch^c5 -c6).*exp(-c7.*weird_constant) ;
    CP=max(0,CP);
    plot3(pitch,tsr,CP);
    grid on;
    hold on;

end


%Graph
figure();
plot(tsr,cp1,'LineWidth',2);hold on;grid on;
plot(lam,cp,'ko');
%plot(lam,Cpp,'ko');
for ii=1:1:3
 txt{ii}=['v_w = ' num2str(vw1(ii)) ' m/s'];
 text(lam(ii),cp(ii),txt{ii},'FontSize',18);
end;
xlabel('lamda','FontSize',14);
xlabel('C_p','FontSize',14);
%% Obtention of the REAL turbine [tipspeedratio(TSR)−Cp] curve
% CP vs tsr curve creation
figure();

vw4=0.1:.1:55;
tsr2= min(20,max(w_ini*(rotor_diameter/2) ./ vw4,0));
k2=(tsr2+c8*angle_pitch).^(-1)-c9/(1+angle_pitch^3);                             %Auxvariable to calculate Cp
cp2=max(0,c1*(c2*k2-c3*angle_pitch-c4*angle_pitch^c5-c6).*exp(-c7*k2));          %calculation of the Cp for←- different tsrs
P2=0.5*rho*cp2*(pi*(rotor_diameter^2)/4).*vw4.^3;
P2a=min(0.5*rho*cp2*(pi*(rotor_diameter^2)/4).*vw4.^3,power);
cp2_new=2.*P2a ./(rho*(pi*(rotor_diameter^2)/4).*vw4.^3);

plot(tsr2,cp2);hold on;grid on;
plot(tsr2,cp2_new);hold on;
plot(lam,cp_new,'ko');
%plot(lam,Cpp,'ko');
for ii=1:1:3
 txt{ii}=['v_w = ' num2str(vw1(ii)) ' m/s'];
 text(lam(ii),cp_new(ii),txt{ii},'FontSize',18);
end;
xlabel('\lambda','FontSize',18);
ylabel('C_p','FontSize',18);


%% MECHANICAL AND ELECTRICAL ANALYSIS OF FIX SPEED WIND TURBINES.
% Parameters
N=transmission_ratio;       %ratio
D=rotor_diameter;           %diameter
poles=pole_pairs;           %pole pairs
Ugrid=[960,850,960,960];      %grid voltage    
f=[50,50,53,47];                %freq
vww=[7,11,14];

% Equivalent model
Zs=Rs+j*Xs;   Zm=Rm*(j*Xm)/(Rm+j*Xm);    Zr=Rr+j*Xr;
Zsm=Zs*Zm/(Zs+Zm);
Zth=Zr+Zsm;


%RESULT MATRICES
Pmm=zeros(4,3);
Tmm=zeros(4,3);
slipp=zeros(4,3);
welec=zeros(4,3);
wsync=zeros(4,3);
Nsync=zeros(4,3);
wgen=zeros(4,3);
Ngen=zeros(4,3);
wturb=zeros(4,3);
Nturb=zeros(4,3);
Igrid=zeros(4,3);
Iabs=zeros(4,3);
Iangle=zeros(4,3);
Sgrid=zeros(4,3);
Pgrid=zeros(4,3);
Qgrid=zeros(4,3);
losses=zeros(4,3);
efficiency=zeros(4,3);
Cpfinal=zeros(4,3);
pitch=zeros(4,3);
lamda=zeros(4,3);
for i=1:4

V=Ugrid(i)/sqrt(3);
Uth= Zm*V/(Zm+Zs);
Uu=abs(Uth);

for j=1:3
angle_pitch=0;   
vw=vww(j);

% Initialization
ss=0;kk=0;err=1;P1=0;
while ((kk<100)&&(err>1e-6))
 kk=kk+1;
 we=2*pi*f(i);  ws=we/poles;  wg=ws*(1-ss);  wt=wg/N;
 %lam= wt*(D/2) / vw;
 %Cpp= 0.0045 * (100 - (lam-10)^2);

 tsr3= min(20,max(wt*(D/2) / vw,0));
 k3=(tsr3+c8*angle_pitch)^(-1)-c9/(1+angle_pitch^3);                             %Auxvariable to calculate Cp
 cp3=max(0,c1*(c2*k3-c3*angle_pitch-c4*angle_pitch^c5-c6)*exp(-c7*k3));          %calculation of the Cp for←- different tsrs
 %{
 Pmech=0.5*rho*cp3*(pi*(D^2)/4)*vw^3;
 if Pmech>power
    Pmech=min(0.5*rho*cp3*(pi*(D^2)/4)*vw^3,power);
    cp3new=2*Pmech/(rho*(pi*(rotor_diameter^2)/4)*vw^3);
    cc1 = 0.73;
    cc2 =151;
    cc3 =0.58; cc4 =0.002; cc5 =2;
    cc6 = 13.2;
    cc7 = 18 ;
    cc8 =-0.02;
    cc9 = -0.003;
    error=cp3-cp3new;
    while error>0.001
        angle_pitch=angle_pitch+0.1;
        kk3=(tsr3+cc8*angle_pitch)^(-1)-cc9/(1+angle_pitch^3);
        cpp3=max(0,cc1*(cc2*kk3-cc3*angle_pitch-cc4*angle_pitch^cc5-cc6)*exp(-cc7*kk3));
        error=cpp3-cp3new;
    end
    cp3=cp3new;
 end
 %}
 
 Pmech=min(0.5*rho*cp3*(pi*(D^2)/4)*vw^3,power);
 err=abs(Pmech-P1);
 P1=Pmech;
 P1ph=-P1/3; 
 a1=1; a2=2*real(Zth)*P1ph-Uu^2; a3=abs(Zth)^2*(P1ph^2);
 sol=roots([a1 0 a2 0 a3]);
 Rt=sol.^2/P1ph;
 s=Rr./(Rr+Rt);
 if (imag(s)~=0)
     s=-abs(s);
 end
 ss=s(2);
end;
kk
Pmm(i,j)=P1;
Tmm(i,j)=P1/wt;
Cpfinal(i,j)=cp3;
pitch(i,j)=angle_pitch;
lamda(i,j)=tsr3;
slipp(i,j)=ss;
welec(i,j)=we;
wsync(i,j)=ws;
Nsync(i,j)=60*ws/ (2*pi);
wgen(i,j)=wg;
Ngen(i,j)=60*wg/ (2*pi);
wturb(i,j)=wt;
Nturb(i,j)=60*wt/ (2*pi);
Igrid(i,j) = V/(Zs+Zm*(Zr+Rt(2))/(Zm+(Zr+Rt(2))));
Iabs(i,j)= sqrt((real(Igrid(i,j))^2+imag(Igrid(i,j))^2));
Iangle(i,j)=rad2deg(atan(imag(Igrid(i,j))/real(Igrid(i,j))));
Sgrid(i,j)= 3*V*conj(Igrid(i,j));
Pgrid(i,j)=-real(Sgrid(i,j));
Qgrid(i,j)=-imag(Sgrid(i,j));
losses(i,j)=real(Sgrid(i,j))+P1;
efficiency(i,j)=Pgrid(i,j)/P1;

end

end