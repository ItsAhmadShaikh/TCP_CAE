%% Compute orientation of the cocentric cylinders in the twisted fiber
clear all; clc; close all;

rf = 1.0; % fiber outer radius
alpha_f = 35; % biad angle at the outermost surface of the fiber
N = 50; % number of concetric cylinders;
% r = linspace(0,rf,N);
r = [0 0.0625, 0.25, 0.5, 0.75, 1 ];
alpha_r = atand((r./rf).*tand(alpha_f));

plot(r,alpha_r,'r-');
hold on;
plot([0 rf],[0,alpha_f],'k-')

%% Properties of the the precursor fiber from Diego R Higueras Ruiz et al (2019) MS thesis
NT = 50;
T0 = 25;
Tf = 120;
T = linspace(T0,Tf, NT);
TK = T+273;
E11 = (-7.86350662691297e-6.*(T.^3) + 0.00231667099158860.*T.^2 -0.229423055654062.*T + 8.18373059906981).*1000; % in MPa
E22 = 600.*ones(size(T)); %in MPa
G12 = (-0.00017.*T.^3 + 0.1085.*T.^2 - 16.56.*T + 845.65); %in MPa
G23 = 500.*ones(size(T)); % in MPa
nu12 = 0.3.*ones(size(T));
alpha11 = ( -0.0273.*( T.^2 -T0.^2) + 1.91.*(T-T0) - 32.1 ).*1e-6; % per deg C
alpha22 = ( 2.62.*(T-T0) + 9.29 ).*1e-6; % per deg C

figure; plot(T,E11)
figure; plot(T,G12)
figure; plot(T,alpha11)
figure; plot(T,alpha22)

thermalProps = [alpha22' alpha11' alpha22' TK']; % to generate a table for abaqus
ElasticPros = [E22' E11' E22' nu12' nu12' nu12' G12' G23' G12' TK'];
