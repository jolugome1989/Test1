%% Creation of a GSCM model
% This code is for generating a wireless scenario with MIMO antennas. The
% MIMO array is going to be an Uniform Linear Array.
% The scenario for this simulation is going to be scatterers surrounding
% the User Equipment.
% Copyright, Jorge Gomez, @USC 2018
%% Clean MATLAB environment
clear;close all;clc
%% Parameters of the simulation
Mt=8; % Number Tx antennas
Mr=8; % Number Rx antennas
fc=2.4e9; % Carrier frequency
c0=3e8; % Speed of light
Nio=20; % Number of interfering objects
R=30; % radius of the scatterers surrounding UE
rho_power=1;
alpha=2; % Path exponent for LOS
SNR_dB=10;
SNR_linear=10.^(SNR_dB/10);
lambda_0=c0/fc; % Wavelength
d_sc=lambda_0/2;
bs=[-100,0]; % Base station position
ue=[0,0]; % User Equipment
pos_clusters=ue;
n_iter=1e3;
io=scatterers_pos(Nio,bs,'Circle',R,ue,'NLOS'); % Function to generate IO in a GSCM scenario
for pos_delta=1:9
[dist_io,t_p] = calc_dist(io,c0,bs,ue,pos_clusters); % Calculate distance and delay
P_io=(dist_io).^(-alpha); % Calcualte Path Loss for each IO using log-distance power law equation
del_pha=exp(-1j*2*pi/lambda_0.*dist_io);
[ste_ma,AOA,AOD,ste_tx,ste_rx]=ULA_model(io,bs,ue,Mt,Mr,d_sc,lambda_0); % Uniform Linear Array Model for eah IO (MtxMrxNio)
C=zeros(1,n_iter);
H_norm=C;
H_n=zeros(Mt,Mr,n_iter);
for kk=1:n_iter
rand_phase=exp(1j.*unifrnd(-pi,pi,[Nio,1]));
com_amp_io=sqrt(P_io).*rand_phase.*del_pha; % Complex amplitude for each IO  (Niox1)
H_system=zeros(Mt,Mr);
    for jj=1:Nio
        H=com_amp_io(jj)*ste_ma(:,:,jj);
        H_system=H_system+H;
    end
H_norm(kk)=norm(H_system,'fro').^2;
H_n(:,:,kk)=H_system;
end
power_sys=1/(Mt*Mr)*mean(H_norm);% Normalization of H matrix is the problem
for iter=1:n_iter
C(iter)=log2(abs(det(eye(Mr)+SNR_linear/Mt*(1/power_sys)*(H_n(:,:,iter)*H_n(:,:,iter)'))));
end
ue(1)=ue(1)+100;
figure(10)
[cdf_plot,stats]=cdfplot(C), hold on
end
%% Verify the position of the BS, UE and IO
figure(1),plot(io(:,1),io(:,2),'*b'),xlabel('xcoord [m]'),ylabel('ycoord [m]'),title('GSCM Scenario'),grid minor, hold on
plot([bs(1) ue(1) 0],[bs(2) ue(2) 0],'^r'),hold on
plot([bs(1) ue(1) 0],[bs(2) ue(2) 0],'--r')