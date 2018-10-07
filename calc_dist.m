%% Function calc_dist:
% Input: scatterers position, Tx and Rx position, number of Tx and Rx
% antennas and center frequency
% Output: travelled distance by each ray, delay of each ray, Tx and Rx steering
% matrix, AoA and AoD of each ray
% For the steering matrix, the distance between antenna elements is
% lambda/2
function [dist_total,t_p,distx_mpc,disrx_mpc,dist_cluster] = calc_dist(pos_scat,c0,pos_tx,pos_rx,cen_cluster)

distx_mpc=sqrt((pos_scat(:,2)-pos_tx(2)).^2+(pos_scat(:,1)-pos_tx(1)).^2);

disrx_mpc=sqrt((pos_rx(2)-pos_scat(:,2)).^2+(pos_rx(1)-pos_scat(:,1)).^2);

dist_total=distx_mpc+disrx_mpc;

t_p=dist_total./c0;

distx_cen_cluster=sqrt((cen_cluster(:,2)-pos_tx(2)).^2+(cen_cluster(:,1)-pos_tx(1)).^2);

disrx_cen_cluster=sqrt((cen_cluster(:,2)-pos_rx(2)).^2+(cen_cluster(:,1)-pos_rx(1)).^2);

dist_cluster=distx_cen_cluster+disrx_cen_cluster;

end