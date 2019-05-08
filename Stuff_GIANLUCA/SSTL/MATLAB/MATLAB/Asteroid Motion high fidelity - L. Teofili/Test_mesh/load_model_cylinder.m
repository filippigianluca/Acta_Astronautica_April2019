% The function takes as input a string that is the position of the file of interest
% e.g. 'data/4179toutatis.tab' and the position of the center of mass and gives as output vertex vectors [r1 r2 r3] and face vectors
% rf_ centered in the baricenter of the asteoid (if no position vector is
% the baricenter will be fixed in the origin).

function [r1,r2,r3,vert,faces,rf_] = load_model_cylinder(filename,R_cm)

       
% fv.vertices = 1*Data.data(1:nvert,:);  % put 1 instead of 1000 if there are problems
% fv.faces = Data.data((nvert+1):end,:);
% 
% vert = fv.vertices;
% faces = fv.faces;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% to load a cilinder

mesh=loadmesh2D('cylinder.msh');
vert = mesh.vert;
l=length(vert);
R_cm = [0,0,-14];
R_cm = repmat(R_cm,l,1);
vert = [vert(:,1)*0.001,vert(:,2)*0.001,0.1*vert(:,3)]+R_cm;
faces = mesh.face;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    
     figure(1)
     fv.facevertexcdata = sqrt(sum(g.^2,2));
%    
     %fv.vertices = 1*Data.data(1:nvert,:)+R_cm; % just to make the plot
     
     fv.vertices = vert;
     fv.faces = mesh.face;
%     
     %patch(fv,'FaceColor','flat')
     patch(fv,'FaceColor','white')
     title('Geographos')%filename)
     axis equal
     xlabel('X [Km]') 
     ylabel('Y [Km]')
     zlabel('Z [Km]')
%     %colorbar
%     %ylabel(colorbar,'acceleration [m/s^2]','FontSize', 16,'FontName', 'Times New Roman')
%     %set(gca, 'FontSize', 16,'FontName', 'Times New Roman')
%     grid on
%     
%     if strcmp('data/massive_point.tab',filename)==1
%         hold on
%         plot3(R_cm(1,1),R_cm(1,2),R_cm(1,3),'o');
%     end


end