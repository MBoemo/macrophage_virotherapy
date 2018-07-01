function [ out ] = fun_p_m_fix(m_coordinate,c_coordinate)

global c_p;
global m_p;
global p_m_max;



out = p_m_max*fun_S(c_coordinate,c_p)*(1 - fun_S(m_coordinate,m_p));

end

