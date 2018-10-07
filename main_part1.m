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
lambda_0=c0/fc; % Wavelength
Nio=100; % Number of interfering objects
bs=[-300*lambda_0,0]; % Base station position
ue=[0,0]; % User Equipment
pos_clusters=ue;
d_sc=lambda_0/2;
rayl_dist=2*(Mt*d_sc)^2/lambda_0;
R=200*lambda_0; % radius of the scatterers surrounding UE
rho_power=1;
alpha=2; % Path exponent for LOS
SNR_dB=3;
SNR_linear=10.^(SNR_dB/10);
io=scatterers_pos(Nio,bs,'Circle',R,pos_clusters,'NLOS'); % Function to generate IO in a GSCM scenario
[dist_io,t_p] = calc_dist(io,c0,ue,bs,pos_clusters); % Calculate distance and delay
P_io=(dist_io).^(-alpha); % Calculate Path Loss for each IO using log-distance power law equation
del_pha=exp(-1j*(2*pi/lambda_0).*dist_io);
for antenna_iter=1:5
d_sc=lambda_0/2^(antenna_iter);
rayl_dist=2*(Mt*d_sc)^2/lambda_0;
[ste_ma,AOA,AOD,ste_tx,ste_rx]=ULA_model(io,bs,ue,Mt,Mr,d_sc,lambda_0); % Uniform Linear Array Model for eah IO (MtxMrxNio)
n_iter=1e4;
C=zeros(1,n_iter);
H_norm=C;
H_n=zeros(Mr,Mt,n_iter);
for kk=1:n_iter
rand_phase=exp(1j.*unifrnd(0,2*pi,[Nio,1]));
com_amp_io=sqrt(P_io).*rand_phase.*del_pha; % Complex amplitude for each IO  (Niox1)
H_system=zeros(Mr,Mt);
    for jj=1:Nio
        H=com_amp_io(jj)*ste_ma(:,:,jj);
        H_system=H_system+H;
    end
H_norm(kk)=norm(H_system,'fro').^2;
H_n(:,:,kk)=H_system;
end
power_sys=1/(Mt*Mr)*mean(H_norm);% Normalization of H matrix is the problem
for iter=1:n_iter
C(iter)=abs(log2(det(eye(Mr)+(SNR_linear/Mt)*(1/power_sys)*(H_n(:,:,iter)*H_n(:,:,iter)'))));
end
figure(10)
[cdf_plot,stats]=cdfplot(C);
xlabel('Rate[bps/Hz]'); ylabel('CDF');hold on
end

%% Verify the position of the BS, UE and IO
figure(1),plot(io(:,1),io(:,2),'*b'),xlabel('x_{coord} [m]'),ylabel('y_{coord} [m]'),title('GSCM Scenario'),grid minor, hold on
plot([bs(1) ue(1)],[bs(2) ue(2)],'^r'),hold on
plot([bs(1) ue(1)],[bs(2) ue(2)],'--r')
% figure(2)
% subplot(2,1,1)
% histogram(AOA,'Normalization','pdf'),grid minor,title('AOA')
% subplot(2,1,2)
% histogram(AOD,'Normalization','pdf'),grid minor,title('AOD')
%stats
%figure(11)
%Ergodic_Capacity_CDF


