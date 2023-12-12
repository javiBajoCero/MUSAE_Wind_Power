clear;
close all;
clc;
%%
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % Dynamic analysis of a fixed −speed wind turbine %
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Wind turbine data (2MW)
% Aerodynamic and mechanic parameters of the wind turbine 
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

%Electrical parameters of the wind turbine (generator)  
power               =2e6;   %[W]
nominal_voltage     =960;   %[Vph-ph]
nominal_current     =1300;  %[A]
connection          ="delta";
pole_pairs          =2;     %[]
rated_speed_50hz    =1500;  %[1/min]
moment_of_inertia   =90;    %[kgm2]

stator_resistance   =0.005; %[ohms]
stator_leakage_inductance   =4e-4;%[H]
rotor_resistance    =0.009;%[ohms]
rotor_leakage_inductance=3e-4;%[H]
magnetizing_inductance=15e-3;%[H]
iron_branch_equivalent_resistance=140;%[ohms]

% Cp values [ Table ]
c1 = 0.44;
c2 =125;
c3 =0; c4 =0; c5 =0;
c6 = 6.94 ;
c7 = 16.5 ;
c8 =0;
c9 = -0.002;

%% simulation and plots
simtime=10;
out=sim("fixed_speed_turbine_simulation.slx");
