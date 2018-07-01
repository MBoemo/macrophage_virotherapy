function [ out ] = fun_d_l(c_coordinate)

global d_l_max;
global c_c;

out = d_l_max*(1 - fun_S(c_coordinate,c_c));


end

