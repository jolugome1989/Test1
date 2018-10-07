N = 1000;
NBS = 2;
NUE = 2;
Ns = 100;
c = 3*10^8;
f = 10^9;
lamda = c/f;
BS = [0; 0];
UE = [120; 120];
rho = 1;
SNR_dB = 3;
SNR = 10^(SNR_dB/10);
d = 0.5*lamda;
alpha = 2;

from = 80;
to = 120;
S = (to-from).*rand(2,Ns) + from;

d1 = sqrt((S(2,:)-UE(2,1)).^2+(S(1,:)-UE(1,1)).^2);
d2 = sqrt((S(2,:)-BS(2,1)).^2+(S(1,:)-BS(1,1)).^2);

DOD = atan2(S(2,:)-UE(2,1),S(1,:)-UE(1,1));
DOA = atan2(S(2,:)-BS(2,1),S(1,:)-BS(1,1));

H = zeros(NBS,NUE,N);
for i = 1:Ns
    for k = 1:NUE
        for m = 1:NBS
            for n = 1:N
                phi = (2*pi)*rand();
                a = rho*exp(1i*phi)/((d1(i)+d2(i))^(alpha/2));
                delay_shift = exp(-1i*2*pi*(d1(i)+d2(i))/lamda);
                DOA_shift = exp(-1i*2*pi*d*(m-1)*sin(DOA(i))/lamda);
                DOD_shift = exp(-1i*2*pi*d*(k-1)*sin(DOD(i))/lamda);
                H(m,k,n) = H(m,k,n) + a*delay_shift*DOA_shift*DOD_shift;
            end
        end
    end
end

H_norm = zeros(1,N);
for n = 1:N
    H_norm(n) = norm(H(:,:,n),'fro').^2;
end

H_avg = mean(H_norm)/(NBS*NUE);

C = zeros(1,N);
for n = 1:N
    C(n) = abs(log2(det(eye(NBS)+(SNR/NUE)*(1/H_avg)*(H(:,:,n)*H(:,:,n)'))));
end

cdfplot(C);