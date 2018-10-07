function[io_pos]=scatterers_pos(Nio,bs,pdf_pos,R,pos_clusters,scenario)
[ncluster,~]=size(pos_clusters);
io_pos=[];
if strcmp(pdf_pos,'Laplacian')
    for j=1:ncluster
    cluster=zeros(Nio,2);
    for i=1:Nio
    cluster(i,:) = mvlaprnd(2,pos_clusters(j,:)',[R^2 0;0 R^2]);
    end
    end
    io_pos=[io_pos;cluster];
elseif strcmp(pdf_pos,'Gaussian') % Gaussian otherwise
    for i=1:ncluster % Create scatterers position using a Jointly Gaussian Distribution
    cluster=mvnrnd(pos_clusters(i,:),[R^2,0;0 R^2],Nio);
    io_pos=[io_pos;cluster];
    end
else
    for i=1:ncluster % Create scatterers position in a circle of radius R
    angles=2*pi.*rand([Nio,1]);
    cluster=[R.*cos(angles),R.*sin(angles)]+pos_clusters(i,:);
    io_pos=[io_pos;cluster];
    end
end
% if strcmp(scenario,'LOS')
% io_pos=[bs;io_pos];
% end
end
