function mesh = loadmesh2D(filename)
switch nargin
    case 0
    filename = 'cylinder.msh'
    otherwise
        
end

fid = fopen(filename,'r');

cur_line = fgetl(fid);
% vertex, face and element index
vi=1;
fi=1;
ei=1;
cur_line = fgetl(fid);
while ~strcmp(cur_line,'$EndElements')
    cur_line = fgetl(fid);
    line = str2num(cur_line);
    if length(line)==4
        mesh.vert(vi,(1:3))=line(2:4);
        vi = vi+1;
    end
    if length(line)==8
        mesh.face(fi,(1:3))=line(6:8);
        fi = fi+1;
    end
end
mesh.nvert = vi-1;
mesh.nfaces = fi-1;
mesh.nelem = ei-1;
end