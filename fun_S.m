function [ out ] = fun_S(c_coordinate,c_constant)

global alpha;

out = (c_coordinate^alpha*(c_constant^alpha + 1))/(c_constant^alpha + c_coordinate^alpha);


end

