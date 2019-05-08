function out = R2fun(TIME_vector, N_images_vector, tt)

% global N_images_vector2R_function
% N_images_vector = N_images_vector2R_function;
% 
% global TIME_vector2R_function
% TIME_vector = TIME_vector2R_function;


out = interp1(TIME_vector, N_images_vector, tt,'linear','extrap');
end