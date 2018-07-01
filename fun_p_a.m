function [ out ] = fun_p_a(l_coordinate,m_coordinate,c_coordinate)

global c_p;
global p_am_max;
global p_al_max;

out = (1 - fun_S(c_coordinate,c_p))*(p_al_max*l_coordinate + p_am_max*m_coordinate);


end

