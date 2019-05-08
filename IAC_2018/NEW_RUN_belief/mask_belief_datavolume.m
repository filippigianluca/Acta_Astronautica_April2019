function res = mask_belief_datavolume(d, u, par)

output = Acta_Astronautica_objconstr1(d, u, par);
 
res =  -(par.fix.nu - output.c);
end