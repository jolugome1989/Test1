function[ste_ma,AOA,AOD,ste_tx,ste_rx]=ULA_model(io,bs,ue,Mt,Mr,d_sc,lambda_0) % Uniform Linear Array Model
[num_io,~]=size(io);
AOD=atan2((io(:,2)-ue(2)),(io(:,1)-ue(1)));
AOA=atan2((io(:,2)-bs(2)),(io(:,1)-bs(1)));


ste_tx=exp(-1j*(2*pi/lambda_0)*d_sc.*(0:Mr-1).*sin(AOD));
ste_rx=exp(-1j*(2*pi/lambda_0)*d_sc.*(0:Mt-1).*sin(AOA));
ste_ma=zeros(Mr,Mt,num_io);
for kk=1:Mr
    for jj=1:Mt
        for ii=1:num_io
            ste_ma(kk,jj,ii)=exp(-1j*(2*pi/lambda_0)*d_sc*(kk-1)*sin(AOD(ii)))*exp(-1j*(2*pi/lambda_0)*d_sc*(jj-1).*sin(AOA(ii)));
        end
    end
end

end