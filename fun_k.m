function [ out ] = fun_k(c_coordinate)

global k_max;
global c_p;

out = k_max*(1 - fun_S(c_coordinate,c_p));


end

