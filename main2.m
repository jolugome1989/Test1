%% Creation of a GSCM model
% This code is for generating a wireless scenario with MIMO antennas. The
% MIMO array is going to be an Uniform Linear Array.
% The scenario for this simulation is going to be scatterers surrounding
% the User Equipment.
% Copyright, Jorge Gomez, @USC 2018
%% Clean MATLAB environment
clear;close all;clc
%% Parameters of the simulation
Mt=4; % Number Tx antennas
Mr=4; % Number Rx antennas
fc=5e9; % Carrier frequency
c0=3e8; % Speed of light
Nio=200; % Number of interfering objects
bs=[0,0]; % Base station position
ue=[50,0]; % User Equipment
pos_clusters=ue;
lambda_0=c0/fc; % Wavelength
d_sc=lambda_0*5;
R=500; % radius of the scatterers surrounding UE
rho_power=1;
alpha=2; % Path exponent for LOS
SNR_dB=10;
SNR_linear=10.^(SNR_dB/10);
C=zeros(1,1000);
H_norm=C;
H_n=zeros(Mt,Mr,1000);
for k=1:1000
io=scatterers_pos(Nio,bs,'Circle',R,pos_clusters,'NLOS'); % Function to generate IO in a GSCM scenario
[dist_io,t_p] = calc_dist(io,c0,bs,ue,pos_clusters); % Calculate distance and delay
P_io=(dist_io).^(-alpha); % Calcualte Path Loss for each IO using log-distance power law equation
del_pha=exp(-1j*2*pi/lambda_0.*dist_io);
[ste_ma,AOA,AOD,ste_tx,ste_rx]=ULA_model(io,bs,ue,Mt,Mr,d_sc,c0,fc); % Uniform Linear Array Model for eah IO (MtxMrxNio)
com_amp_io=sqrt(P_io).*del_pha; % Complex amplitude for each IO  (Niox1)
H_system=zeros(Mt,Mr);
    for j=1:Nio
        H=com_amp_io(j)*ste_ma(:,:,j);
        H_system=H_system+H;
    end
H_norm(k)=norm(H_system,'fro').^2;
H_n(:,:,k)=H_system;
end
power_sys=1/(Mt*Mr)*mean(H_norm);% Normalization of H matrix is the problem
H_nn=H_n/sqrt(power_sys);
for iter=1:1000
C(iter)=log2(abs(det(eye(Mr)+SNR_linear/Mt*(H_nn(:,:,iter)*H_nn(:,:,iter).'))));
end

%% Verify the position of the BS, UE and IO
% figure(1),plot(io(:,1),io(:,2),'*b'),xlabel('xcoord [m]'),ylabel('ycoord [m]'),title('GSCM Scenario'),grid minor, hold on
% plot([bs(1) ue(1)],[bs(2) ue(2)],'^r'),hold on
% plot([bs(1) ue(1)],[bs(2) ue(2)],'--r')
% figure(2)
% subplot(2,1,1)
% histogram(AOA,'Normalization','pdf'),grid minor,title('AOA')
% subplot(2,1,2)
% histogram(AOD,'Normalization','pdf'),grid minor,title('AOD')
figure(10)
[cdf_plot,stats]=cdfplot(C);
stats
figure(11)
Ergodic_Capacity_CDF