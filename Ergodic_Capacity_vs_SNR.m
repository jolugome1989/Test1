% Ergodic_Capacity_vs_SNR.m
clear, close all
SNR_dB=0:5:20; SNR_linear=10.^(SNR_dB/10);
N_iter=1000; sq2 = sqrt(0.5);
for Icase=1:6
if Icase==1, nT=1; nR=1; % 1x1
elseif Icase==2, nT=1; nR=2; % 1x2
elseif Icase==3, nT=2; nR=1; % 2x1
elseif Icase==4, nT=2; nR=2; % 2x2
elseif Icase==5, nT=4; nR=4; % 4x4
else
    nT=8; nR=8; % 8x8
end
n=min(nT,nR); I = eye(n);
C(Icase,:) = zeros(1,length(SNR_dB));
for iter=1:N_iter
H = sq2*(randn(nR,nT)+1j*randn(nR,nT));
if nR>=nT
    HH = H'*H; 
else 
    HH = H*H'; 
end
for i=1:length(SNR_dB) % Random channel generation
C(Icase,i) = C(Icase,i)+log2(real(det(I+SNR_linear(i)/nT*HH)));
end
end
end
C = C/N_iter;
plot(SNR_dB,C(1,:),'b-o', SNR_dB,C(2,:),'b-', SNR_dB,C(3,:),'b-s');
hold on, plot(SNR_dB,C(4,:),'b->', SNR_dB,C(5,:),'b-^',SNR_dB,C(6,:),'b-<','linewidth',2);grid minor
xlabel('SNR[dB]'); ylabel('bps/Hz');