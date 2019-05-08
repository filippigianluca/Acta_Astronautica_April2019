function [F] = Fsstl2(d);

    F(1) = sstl2c(d);
    F(2) = 1-sstl2p(d);

return