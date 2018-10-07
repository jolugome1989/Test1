% Rui Wang
% Date: May 18th, 2013

% Function: Computes the optimal phase vector for multitone signals with 
% low crest factor

% update the function to select the bandwidth and upper and lower
% freqeucnies

% Reference: Mathias Friese, Multitone Signals with Low Crest Factor, IEEE
% Trans on Communications, Vol. 45, No.10, Oct. 1997

%{
Revision History
Programmer                      Change                     Date
Rui Wang                       original                May 18th,2013
Rui Wang                  change parameters for Usrp   Nov 24th,2015
Rui Wang                add comments for f1 = 0        May 5th, 2016
Rui Wang         modify the script to ensure more ifft points and better
                 accuracy for the time-domain signal   Nov 8th, 2017
Jorge Gomez      Change the frequency sampling to 80 MS/s  Dec 10th, 2017
                 plot the real and imaginary parts      
%}
% Copyright Rui Wang 2013

clear;clc;
close all;

%Nf = 384; % number of frequency tones
N = 4000; %N = 80; % number of IFFT points
fs = 100e6;  % sampling frequency Hz
T0 = 1e-6; % OFDM symbol duration
df = 1/T0; % subcarrier spacing in Hz
f1 = 0*df; fix1 = floor(f1/df);% low frequency, index
f2 = 40e6; fix2 = floor(f2/df);% high frequency, index

if f1 == 0
    % case DC component is included, f1 = 0
    Nf = 2*(fix2 - fix1 + 1) -1; % odd number
    %Nf = 2*fix2 + 1;
    ftone = zeros(1,Nf); % frequency tones
    ftone(1:floor(Nf/2)) = (-fix2*df):df:-df; % negative frequency tones
    ftone(floor(Nf/2)+1) = 0; % DC 
    ftone(ceil(Nf/2)+1:end) = df:df:fix2*df; % positive frequency tones
    fix_eff = 1:Nf; % useful subcarrier index
else
    Nf_tot = 2*(fix2 - fix1 + 1); % even number
    ftone = zeros(1,Nf_tot);
    ftone(1:Nf_tot/2) = (-fix2*df):df:(-fix1*df);
    ftone(Nf_tot/2+1:Nf_tot) = -ftone(Nf_tot/2:-1:1);
    %ftone(Nf/2+1:Nf) = -ftone(1:Nf/2); % old incorrect
    Nf = 2*fix2+1;
    fix_eff = [1:Nf_tot/2,Nf-Nf_tot/2+1:Nf];
end

% normal time index
t = 0:1/fs:(fs/df-1)/fs; 
% densified time index
t_D = linspace(0,T0,N+1); t_D = t_D(1:end-1); 


% n = zeros(1,Nf);
% if mod(Nf,2) 
%     % when Nf is odd, or f1 = 0
%     %n(1:ceil(Nf/2)) = (-floor(f2/f_res)):1:(-floor(f1/f_res));
%     %n(Nf/2+1:Nf) = -n(1:Nf/2); % old
%     %n(ceil(Nf/2)+1:Nf) = -n(ceil(Nf/2)-1:-1:1);
%     n = -(Nf-1)/2 : (Nf-1)/2;
% else
%     % when Nf is even, or f1 > 0
%     n(1:Nf/2) = -fix2:-fix1;
%     n(Nf/2+1:Nf) = -n(Nf/2:-1:1);    
% end

n = -(Nf-1)/2 : (Nf-1)/2;
Phi = pi*n.^2/Nf;
w_res = 2*pi*df;


ix = 0; % iteration index
max_ix = 10000;
CF = 5;
S = sqrt(Nf);
CF_store = zeros(1,max_ix);
key_reconstruct = 0;
key_liftthreshold = 1;
key_plotm_t = 0;
while (ix<=max_ix) && (S <sqrt(Nf)*CF)
    ix = ix+1;
% Calculate time-domain signal     
%     m_t = zeros(size(t));
%     for kk = 0:Nf-1
%         m_t = m_t+exp(1j*(2*pi*ftone(kk+1)*t+Phi(kk+1)));
%     end

    M = exp(1j*Phi);
    if length(fix_eff) ~= Nf
        M(setdiff(1:Nf,fix_eff)) = zeros(1,Nf-length(fix_eff));
    end
  
    M_sft = ifftshift(M);
    % pad zeros in the middle frequency, i.e. after the highest frequency
    M_pad = [M_sft(1:ceil(Nf/2)),zeros(1,N-Nf),M_sft(ceil(Nf/2)+1:end)]; 
    m_t = N*ifft(M_pad);
    if key_plotm_t ==1
        figure(10)
        plot(abs(m_t));
        figure(11)
        plot(abs(fft(m_t)));
    end
    CF_store(ix) = max(abs(m_t))/sqrt(mean(abs(m_t).^2));
    CF = CF_store(ix);
    
% Calculate error signal
    e_t = (abs(m_t)>=S).*(abs(m_t)-S).*exp(1j*angle(m_t));
% Fourier series expansion (OBSOLETE)
%     e_k = zeros(1,Nf);
%     n = 0:N-1;
%     for ii=1:Nf
%         e_k(ii) = 1/N*sum(e_t.*exp(-1j*2*pi*ftone(ii)*n/fs));
%     end
    
    % Find frequency coefficients of e_k
    e_k_pad = 1/N*fft(e_t); 
    e_k_sft = [e_k_pad(1:ceil(Nf/2)),e_k_pad(ceil(Nf/2)+1+N-Nf:end)];
    e_k = fftshift(e_k_sft);
    
    if key_reconstruct == 1
        e_t_recon = zeros(size(t));
        for kk = 1:Nf
            e_t_recon = e_t_recon + e_k(kk)*exp(1j*2*pi*ftone(kk)*t);
        end
        figure(99)
        plot(t,abs(e_t),'r');
        hold on;
        plot(t,abs(e_t_recon),'b');
        legend('original error waveform','renconstruction');
        grid on;
        err_recon = mean(abs(e_t_recon-e_t).^2);
        err_pratio = err_recon/mean(abs(e_t).^2);
    else
%         disp('no reconstruction comparison')
    end
    
% subtraction of error signal in the frequency domain
    for kk = 1:Nf
        Phi(kk) = angle( exp(1j*Phi(kk))-e_k(kk) );
    end

% update the threshold
    if key_liftthreshold == 1
     S = S+0.001;
    else 
        S=S+0;
    end
end    

% Final waveform sampled at fs 
m_t_fs = m_t(1:length(t_D)/length(t):end);
m_t_fs_norm=m_t_fs/max(abs(m_t));

figure(10)
subplot(2,1,1)
plot(t_D*1e6,real(m_t)/max(abs(m_t)));
hold on;
plot(t*1e6,real(m_t_fs_norm),'*','Linewidth',2);grid minor
title('Basis waveform real m(t)'); 
xlabel('t/ \mu s')
legend(sprintf('%dx Upsampled',length(t_D)/length(t)),'Normal Rate');
axis([0,max(t_D*1e6),-1.5,1.5]);
hold off;
subplot(2,1,2)
plot(t_D*1e6,imag(m_t)/max(abs(m_t)));
hold on;
plot(t*1e6,imag(m_t_fs_norm),'*','Linewidth',2);grid minor
title('Basis waveform imaginary m(t)'); 
xlabel('t/ \mu s')
legend(sprintf('%dx Upsampled',length(t_D)/length(t)),'Normal Rate');
axis([0,max(t_D*1e6),-1.5,1.5]);
hold off;

figure(1)
subplot(2,1,1)
plot(t_D*1e6,abs(m_t)/max(abs(m_t)));
hold on;
plot(t*1e6,abs(m_t_fs_norm),'*','Linewidth',2);grid minor
title('Basis waveform |m(t)|'); 
xlabel('t/ \mu s')
legend(sprintf('%dx Upsampled',length(t_D)/length(t)),'Normal Rate');
axis([0,max(t_D*1e6),0,1.5]);
hold off;

f = -fs/2:df:(fs/2-df);
subplot(2,1,2);
plot(f,abs(ifftshift(fft(m_t_fs_norm))),'Linewidth',2);grid minor
xlabel('f/Hz');ylabel(' / dB');
title('Mag of FFT of the basis waveform');

% figure(3)
% f_long = -fs/2:df/2:fs/2-df/2;
% f_ix = 1:2:length(f_long);
% subplot(2,1,1)
% plot(f_long,10*log10(fftshift(abs(fft(m_t,2*N)))));
% xlabel('f/Hz'); ylabel('|mag| / dB');
% title('FFT of 2N samples, zero-padding');
% subplot(2,1,2)
% m_mag = fftshift(abs(fft(m_t,2*N)));
% plot(f,10*log10(m_mag(f_ix)));
% xlabel('f/Hz'); ylabel('|mag| / dB');
% 
% figure(2)
% plot(CF_store(1:ix),'Linewidth',2);grid minor
% xlabel('Iteration number'); ylabel('Crest Factor');
% 
% figure(3) 
% m_rep4 = repmat(m_t,1,4);
% pwelch(m_rep4,128,64,length(m_rep4),fs,'twosided');
% title('4 times the length of original sequence');
% 
% figure(4) 
% m_rep4_wcp = [m_rep4(end-199:end),m_rep4];
% length(m_rep4_wcp)
% pwelch(m_rep4_wcp,128,120,length(m_rep4_wcp),fs,'twosided');
% title('')
save('1MIMO_100MS_80MHz_100pts.mat','m_t_fs','fs','m_t_fs_norm')
M_f_norm=fftshift(fft(m_t_fs_norm));
M_f=fftshift(fft(m_t_fs));
