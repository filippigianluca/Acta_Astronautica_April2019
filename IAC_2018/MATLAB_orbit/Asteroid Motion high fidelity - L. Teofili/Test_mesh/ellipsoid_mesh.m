function [r1,r2,r3,vert,faces,rf_,nfaces] = ellipsoid_mesh(radius1,radius2,radius3,R_cm)

a_I = radius1;
b_I = radius2;
c_I = radius3;

fd=@(p) p(:,1).^2/1+p(:,2).^2/1+p(:,3).^2/1-1;
%[p,t]=distmeshsurface(fd,@huniform,0.2,[-2.1,-1.1,-1.6; 2.1,1.1,1.6]);
%fd=@(p) p(:,1).^2/40^2+p(:,2).^2/1+p(:,3).^2/1^2-1;
[p,t]=distmeshsurface(fd,@huniform,0.1,[-2.1,-1.1,-1.6; 2.1,1.1,1.6]);
%[p,t]=distmeshsurface(fd,@huniform,0.55,[-2.1,-1.1,-1.6; 2.1,1.1,1.6]); %0.55 mesh lassa

p(:,1) = a_I*p(:,1);
p(:,2) = b_I*p(:,2);
p(:,3) = c_I*p(:,3);

fv.vertices = p;

fv.faces = t;

nvert = length(fv.vertices);

if nargin<4
    R_cm = [0,0,0];
end

R_cm = repmat(R_cm,nvert,1);

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

fv.vertices = p + R_cm;

figure(1)
    fv.facevertexcdata = sqrt(sum(g.^2,2));
    patch(fv,'FaceColor','flat')
    %title(filename)
    axis equal
    %colorbar
    %ylabel(colorbar,'acceleration [m/s^2]','FontSize', 16,'FontName', 'Times New Roman')
    %set(gca, 'FontSize', 16,'FontName', 'Times New Roman')
    grid on
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    
end