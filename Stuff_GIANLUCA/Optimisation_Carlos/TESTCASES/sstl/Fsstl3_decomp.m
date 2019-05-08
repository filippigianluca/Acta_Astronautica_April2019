function [F] = Fsstl3_decomp(d);

    F(1) = sstl2c(d);
    F(2) = d(2);
    F(3) = 1-sstl2p_decomp(d);

return