clear all; close all; clc;
%% Machine parameters
Vnom=960;                   % Nominal voltage
Rs=0.005;                   % Stator resistance
Xs=2*pi*50*4e-4;            % Leakage stator inductance (impedance)
Rr=0.009;                   % Rotor resistance
Xr=2*pi*50*3E-4 ;           % Leakage rotor inductance (impedance)
Xm=2*pi*50*15E-3;           % Magnetizing branch inductance(impedance)
Rm=140;                     % Magnetizing branch resistance
V=Vnom/sqrt(3);             % Nominal voltage
pols=2;                     % Pole pairs
ws=2*pi*50/pols;            % Synchronous speed

%% Equivalent scheme
%         Xs      Rs              Xr     Rr/s
% ------|||||--/\/\/\-----------|||||--/\/\/\----.
% +         --->         |  |       --->          |
%            Is          X  /        Ir           |
% V                   Xm X  \                     |
%                        X  / Rm                  |
% -                      |  \                     |
% ------------------------------------------------·
s=-1:.0001:1;                                       %create slip axis swipe
Zm=(Rm*1j*Xm)/(Rm+1j*Xm);                           %magnetizing branch equivalent impedance
paralel=(((Rr./s+1j*Xr)*Zm)./((Rr./s+Xr*1j)+Zm));   %rotor + magnetizing branch equivalent impedance
imp=paralel+Rs+Xs*1j;                               %stator+rotor+magn branch equivalent impedance

Is=V./imp;                                          %stator current
Vr=Is.*paralel;                                     %middle voltage
Ir=Vr./(Rr./s+1j*Xr);                               %rotor current

Telec1=3*abs(Ir.^2)*Rr./(s*ws);                     %Electrical Torque=Power/(s)ws


subplot(3,1,1);
plot(s,Telec1,'LineWidth',2);grid on;
xlabel('s slip [-]','FontSize',14);
ylabel('T_{electric} fast shaft [Nm]','FontSize',14);
legend('Machine Torque')


subplot(3,1,2);
plot((1-s)*ws*30/pi,Telec1,'LineWidth',2);grid on;
xlabel('w_2 fast shaft [rpm]','FontSize',14);
ylabel('T_{elec} fast shaft [Nm]','FontSize',14);
legend('Machine Torque')

%% Turbine parameters
% Cp values [Cp values from the rable]

c1=0.44;
c2=125;
c3=0;c4=0;c5=0;
c6=6.94;
c7=16.5;
c8=0;
c9=-0.002;

% WT parameters
R_turbina=76/2;               % WT radius
A=pi*R_turbina^2;           % Area
rho=1.225;                  % Air density
angle_pitch=0;              % Pitch angle
n_multiplicador = 80;       % Transmission ratio (gearbox)

%% Curve [Speed - turbine torque]
% Calculation for different shaft speed and for different wind speeds

w2=0:1:314;                     % Vector of machine speeds

vw_array=[7,11,14];

for ii=1:1:3                                                                            % Creation of the 'for' loop. It will execute it 6 times, as defined (1:1:6)
    vw=vw_array(ii);                                                                          % Wind speed (for simplicity, we calculate it based on the iteration number)
    clear tsr;                                                                          % Clear the tip speed ratio for the calculation of the current iteration
    tsr= w2*R_turbina/(n_multiplicador*vw);                                             % Calculation of the TSR for the current wind speeds [vw] and the linear blade speed [w2*R_turbina/(n_multiplicador)]
    k1=(tsr+c8*angle_pitch).^(-1)-c9/(1+angle_pitch^3);                                 % aux variable for the cp calculation
    cp=max(0,c1*(c2*k1-c3*angle_pitch-c4*angle_pitch^c5-c6).*exp(-c7*k1));              % Calculation of the Cp constant for the current iteration
    T_turbina(ii,:) = (1/n_multiplicador)*0.5*rho*A*cp*vw^3.*(w2/n_multiplicador).^-1;  % Turbine torque (fast shaft)
    txt{ii}=['Wind=' num2str(vw) ' m/s'];                                               % Creation of a 'txt' vector for the legend
end

%Grafic [Velocitat eix - Parell turbina]
subplot(3,1,3);
plot(w2*30/pi,T_turbina,'LineWidth',2);grid on; % plots the turbine to the speed of the machine
xlabel('\omega_2 fast shaft [rpm]','FontSize',14); % x axis label
ylabel('T_2 fast shaft [Nm]','FontSize',14); % y axis label
legend(txt); % legend
title('Speed - Turbine torque') % title


%% Graphic code
figure();
plot(w2*30/pi,T_turbina,'LineWidth',2);grid on;
hold on
plot((1-s)*ws*30/pi,-Telec1,'LineWidth',2);grid on;
xlabel('\omega_2 fast shaft [rpm]','FontSize',14);
ylabel('T_2 fast shaft [Nm]','FontSize',14);
legend(txt); 
title('Speed - Turbine torque') % title