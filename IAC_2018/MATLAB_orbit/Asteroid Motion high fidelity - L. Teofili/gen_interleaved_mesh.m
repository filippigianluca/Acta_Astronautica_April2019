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

function success = gen_interleaved_mesh(filename, output_name)
       
% if nargin<1
%     filename = 'data/4179toutatis.tab'; %1620geographos.tab';
%     filename = 'data/1998ky26.tab';
% end

Data = importdata(filename);

nvert=0;

while Data.textdata{nvert+1}=='v'
    nvert=nvert+1;
end

verts = 1*Data.data(1:nvert,:);  % put 1 instead of 1000 if there are problems on the asteroid dimention
faces = Data.data((nvert+1):end,:);

norms = cell(size(verts,1),1);
areas = norms;

% Computes face normals and assignes them to each vertex of the respective triangle
for index = 1:size(faces,1)
    A = verts(faces(index,1),:);
    B = verts(faces(index,2),:);
    C = verts(faces(index,3),:);
    
    dir = cross(A-B,A-C);
    
    if (dot(dir,A) < 0) && (dot(dir,B) < 0) && (dot(dir,C) < 0)
        dir = -dir;
    end
    
    area = 0.5*norm(dir);
    dir = dir/norm(dir);
    
    norms{faces(index,1)}{length(norms{faces(index,1)})+1} = dir;
    norms{faces(index,2)}{length(norms{faces(index,2)})+1} = dir;
    norms{faces(index,3)}{length(norms{faces(index,3)})+1} = dir;
    
    areas{faces(index,1)}{length(areas{faces(index,1)})+1} = area;
    areas{faces(index,2)}{length(areas{faces(index,2)})+1} = area;
    areas{faces(index,3)}{length(areas{faces(index,3)})+1} = area;
end

normals = zeros(size(norms,1), 3);
for index = 1:size(norms,1)
    for index2 = 1:size(norms{index},2)
        normals(index,:) = normals(index,:)+areas{index}{index2}*norms{index}{index2};
    end
    normals(index,:) = normals(index,:)/norm(normals(index,:));
end

interleaved_mesh = [verts normals];
faces = faces-1;

fileID = fopen([output_name, '.tab'],'w');
fprintf(fileID,'iv    %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f\n', interleaved_mesh.');        
fprintf(fileID, '\n');
fprintf(fileID,'f    %12i %12i %12i\n', faces.');
fclose(fileID);

disp('blabla')

success = true;