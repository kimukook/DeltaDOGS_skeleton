function [dis] = norm2_dis(x, X)
% calculate the distance from x to each points in X
% x -> n by 1
% X -> n by m
% 
% Author:   Muhan Zhao
% Date  :   May. 20, 2019
dis = sqrt(sum((X - x).^2, 1));
end