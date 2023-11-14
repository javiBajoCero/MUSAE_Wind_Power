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


subplot(2,1,1);
plot(s,Telec1,'LineWidth',2);grid on;
xlabel('s slip [-]','FontSize',14);
ylabel('T_{electric} fast shaft [Nm]','FontSize',14);
legend('Machine Torque')


subplot(2,1,2);
plot((1-s)*ws*30/pi,Telec1,'LineWidth',2);grid on;
xlabel('w_2 fast shaft [rpm]','FontSize',14);
ylabel('T_{elec} fast shaft [Nm]','FontSize',14);
legend('Machine Torque')