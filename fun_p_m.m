function [ out ] = fun_p_m(n_coordinate,c_coordinate)

global c_p;
global n_0;
global p_m_max;
global m_p;

%out = p_m_max*fun_S(c_coordinate,c_p)*min(n_coordinate/n_0,1);

out = p_m_max*fun_S(c_coordinate,c_p)*(1-fun_S(m_coordinate,m_p));


end

