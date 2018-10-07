N_scatter = 100;
N_Ue = 2;
N_Bs = 2;
x_s = zeros(1,N_scatter);
y_s = zeros(1,N_scatter);
angle_depart = zeros(1,N_scatter);
angle_arrive = zeros(1,N_scatter);
f = 2.4*10^9;
c = 3*10^8;
lambda = c/f;
vector = zeros(1,N_scatter);
Distance = 300*lambda;
Bs_x = -Distance;
Bs_y = 0;
R_ns= 200*lambda;%radius
theta = linspace(0,2*pi,200);%you can increase this if this isn't enough yet
x=R_ns*cos(theta);
y=R_ns*sin(theta);
d_S_Bs = zeros(1,N_scatter);


delta_d_Bs = lambda/2;
SNR_Db = 3;
alpha = 2;
SNR = 10.^(0.1*SNR_Db);
N_realization = 1000;
H = zeros(N_Ue,N_Bs,N_realization);
norm_H = zeros(1,N_realization);
C_GSCM = zeros(1,N_realization);
rho = 1;
 
    for Ns = 1:N_scatter
        pos_x = randi(length(x));
        x_s(Ns) = x(pos_x);
        y_s(Ns) = y(pos_x);
        if y_s(Ns)<0
            angle_depart(Ns) = -(180+atan2d(y_s(Ns),x_s(Ns)));
        else 
            angle_depart(Ns) = 180-(atan2d(y_s(Ns),x_s(Ns)));
        end
            angle_arrive(Ns) = atan2d(y_s(Ns)-0,x_s(Ns)-Distance);
    
        d_S_Bs(Ns) = sqrt((y_s(Ns)-Bs_y)^2+(x_s(Ns)-Bs_x)^2);
        d_tot(Ns) = d_S_Bs(Ns)+ R_ns;
    end
for t = 1:N_realization
    for k = 1:N_Ue
        for m = 1:N_Bs
             for Ns = 1:N_scatter
                 H(k,m,t) = H(k,m,t)+ (rho*exp(1i*rand(1)*2*pi)/d_tot(Ns)^(alpha/2)*exp((-1i*2*pi/lambda)*d_tot(Ns))*exp(-1i*(m-1)*cosd(angle_arrive(Ns)))*exp(-1i*(k-1)*cosd(angle_depart(Ns))));
             end
        end
    end
   norm_H(t) = norm(H(:,:,t),'fro').^2;
end
avg_power = mean(norm_H)/(N_Ue*N_Bs);
Ident_mat = eye(N_Bs);
for t = 1:N_realization
    C_GSCM(t) = abs(log2(det(Ident_mat+(SNR/N_Ue)*(1/avg_power)*(H(:,:,t)*H(:,:,t)'))));
end
[H,STATS] = cdfplot(C_GSCM)
%ksdensity(C_GSCM);
%hold on;
%histogram(C_GSCM,'normalization','pdf')