function out = R1fun(TIME_vector, N_images_vector, tt)

% global N_images_vector2R_function
% N_images_vector2R_function = N_images_vector;
% 
% global TIME_vector2R_function
% TIME_vector2R_function = TIME_vector;


out = interp1(TIME_vector, N_images_vector, tt,'linear','extrap');

out=out/2;
end