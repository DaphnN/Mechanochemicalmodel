function J = evaluate_jacobian_num(F,w0,h)
%evaluate_jacobian_num numerically evaluates the Jacobian matrix of the
%provided function F at specified point w0. The derivatives are
%approximated by finite difference scheme
%   Input:
%       F - input function
%       w0 - point, where Jacobian will be evaluated
%   Output:
%       J - Jacobian matrix
% Author: Daphne Nesenberend, Alexey Kazarnikov
% Date: April 2025

M = length(w0);
J = zeros(M,M);
F0 = F(w0);

for j=1:M
    w0h = w0;    
    w0h(j) = w0h(j) + h;
    dF = (F(w0h) - F0)/h;
    J(:,j) = dF;           
end

