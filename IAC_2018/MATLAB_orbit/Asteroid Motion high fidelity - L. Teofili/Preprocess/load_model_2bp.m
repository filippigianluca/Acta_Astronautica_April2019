function [r1,r2,r3,vert,faces,rf_] = load_model_2bp(filename,R_cm)
%load_model_2bp Load asteroid models from database data.
%  [r1,r2,r3,vert,faces,rf_] = load_model_2bp(filename,R_cm) takes as input 
%  a string that is path of the asteroid file of interest
%  e.g. 'data/4179toutatis.tab' and the position of the center of mass and 
%  gives as output vertex vectors [r1 r2 r3], the mesh matrices vert and faces
%  and face vectors rf_ centered in the baricenter position of the asteroid 
%  (if no position vector is specified for the asteroid barycentre it will 
%  be fixed at the origin).
%  Lorenzo Teofili, University of Rome La Sapienza 
%  Nicolas Thiry, University of Strathclyde, 5/11/2015

       
if nargin<1
    filename = 'data/4179toutatis.tab'; %1620geographos.tab';
    filename = 'data/1998ky26.tab';
end

Data = importdata(filename);

nvert=0;

while Data.textdata{nvert+1}=='v'
    nvert=nvert+1;
end

if nargin<2
    R_cm = [0,0,0];
end

R_cm = repmat(R_cm,nvert,1);

fv.vertices = 1*Data.data(1:nvert,:);  % put 1 instead of 1000 if there are problems on the asteroid dimention
fv.faces = Data.data((nvert+1):end,:);

vert = fv.vertices;
faces = fv.faces;

nfaces = length(faces);

rf_ = zeros(nfaces,3);
r1 = rf_;
r2 = r1;
r3 = r1;
g = r1;

for i=1:nfaces
    % vector that identifies a surface on the asteroid 
    rf_(i,:) = (vert(faces(i,1),:) + vert(faces(i,2),:) + vert(faces(i,3),:))/3;
    r1(i,:) = vert(faces(i,1),:);
    r2(i,:) = vert(faces(i,2),:);
    r3(i,:) = vert(faces(i,3),:);
end
    
    
     fv.facevertexcdata = sqrt(sum(g.^2,2));
%    
     fv.vertices = 1*Data.data(1:nvert,:)+R_cm; % just to make the plot
     
     
%     
%       patch(fv,'FaceColor','white')
%       h=patch(fv,'FaceColor',[0.5,0.5,0.5],'Edgecolor','none')
%  h.AmbientStrength =0.25;
%  h.DiffuseStrength = 1;
%  h.SpecularStrength = 0.25;
%  h.SpecularExponent = 25;
%  h.BackFaceLighting = 'lit';
%  lightangle(-90,0)
%  axis equal
%  h=patch(fv,'FaceColor',[0.5,0.5,0.5],'Edgecolor','red','LineWidth',0.01)
     %title('Geographos')%filename)
%      axis equal
%      xlabel('X [Km]') 
%      ylabel('Y [Km]')
%      zlabel('Z [Km]')
% %     %colorbar
%     %ylabel(colorbar,'acceleration [m/s^2]','FontSize', 16,'FontName', 'Times New Roman')
%     %set(gca, 'FontSize', 16,'FontName', 'Times New Roman')
%     grid on
%     
%     if strcmp('data/massive_point.tab',filename)==1
%         hold on
%         plot3(R_cm(1,1),R_cm(1,2),R_cm(1,3),'o');
%     end


end