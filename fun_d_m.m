function [ out ] = fun_d_m(c_coordinate)

global d_m_max;
global c_c;

out = d_m_max*(1 - fun_S(c_coordinate,c_c));

end

