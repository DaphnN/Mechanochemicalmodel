% globalarclenght computes the arclenght of the curve
%   Input:
%       u - height vector
%       N - amount of discretization points
%       h - spatial step
%   Output:
%       Clenght - arclenght of curve
% Author: Daphne Nesenberend, Alexey Kazarnikov
% Date: April 2025

function Clength=globalarclength(u,N,h)

C1=0;
for i = 1:N-1
    C = C1+sqrt(h^2+(u(i)-u(i+1))^2);
    C1=C;
end
Clength=C+sqrt(h^2+(u(1)-u(N))^2);
