
clear;clc;close all;

% WT parameters
R_turbine=45;   %WT radius

% % constant speed (pitch angle=0)
% c1=0.44;
% c2=152;
% c3=0;c4=0;c5=0;
% c6=6.94;
% c7=16.5;
% c8=0;
% c9=-0.002;

% variable speed (variable pitch angle)
c1=0.73;
c2=151;
c3=0.58;c4=0.002;c5=2.14;
c6=13.2;
c7=18.4;
c8=-0.02;
c9=-0.003;


%calculate CP
lambda=1:0.5:15;

for angle_pitch=0:1:29
    pitch=ones(1,29)*angle_pitch;
    weird_constant=1./(lambda+c8.*angle_pitch)-(c9./(1+angle_pitch^3));
    CP=c1.*(c2.*weird_constant -c3.*angle_pitch -c4.*angle_pitch^c5 -c6).*exp(-c7.*weird_constant) ;
    CP=max(0,CP);
    plot3(pitch,lambda,CP);
    grid on;
    hold on;

end