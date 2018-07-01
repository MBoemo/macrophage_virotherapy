function [ out ] = fun_d_phi(c_coordinate)

global c_phi;

out = (1 - fun_S(c_coordinate,c_phi));

end

