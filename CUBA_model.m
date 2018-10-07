function[ste_ma,AOA,AOD,ste_tx,ste_rx]=CUBA_model(io,bs,ue,Mt,Mr,lambda_0,R) % Uniform Linear Array Model

[num_io,~]=size(io);
AOD=atan2((io(:,2)-ue(2)),(io(:,1)-ue(1)));
AOA=atan2((io(:,2)-bs(2)),(io(:,1)-bs(1)));


ste_tx=exp(-1j*2*pi*R/lambda_0.*cos(AOD-2*pi.*(0:(Mr-1))./Mr));
ste_rx=exp(-1j*2*pi*R/lambda_0.*cos(AOA-2*pi.*(0:(Mt-1))./Mt));
ste_ma=zeros(Mr,Mt,num_io);
for kk=1:Mr
    for jj=1:Mt
        for ii=1:num_io
            ste_ma(kk,jj,ii)=exp(-1j*2*pi*R/lambda_0.*cos(AOD(ii)-2*pi.*(kk-1)./Mr))*exp(-1j*2*pi*R/lambda_0.*cos(AOA(ii)-2*pi.*(jj-1)./Mt));
        end
    end
end

end