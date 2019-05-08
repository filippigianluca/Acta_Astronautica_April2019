function bpa = find_bpa(pos_FEc_minmax, bpa_structure)

bpa = 1;
for jj = 1:length(pos_FEc_minmax)
    bpa = bpa*bpa_structure{jj}(pos_FEc_minmax(jj));
end

return