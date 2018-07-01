function [ out ] = fun_d_c(m_coordinate,l_coordinate,n_coordinate,c_coordinate)

global d_cl_max; 
global d_cm_max;
global d_cp_max;

out = c_coordinate*(d_cl_max*l_coordinate + d_cm_max*m_coordinate) + d_cp_max*fun_p_m_fix(m_coordinate,c_coordinate)*m_coordinate;

end

