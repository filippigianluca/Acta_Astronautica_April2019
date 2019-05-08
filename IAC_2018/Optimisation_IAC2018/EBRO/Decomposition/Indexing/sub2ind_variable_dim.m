function ind = sub2ind_variable_dim(matrix_pos, size_matrix)


matrix_pos2ind = '';
for j = 1:length(matrix_pos)
    matrix_pos2ind = strcat(matrix_pos2ind, ',',num2str(matrix_pos(j)));
end


if length(size_matrix) > 1
    ind  = eval(strcat('sub2ind([',num2str(size_matrix),']', matrix_pos2ind, ');'));
else
    ind = matrix_pos;
end
return