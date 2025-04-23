% f computes the curvature dependent morphogen production
%   Input:
%       K - curvature vector 
%       eta - parameter value eta
%       delta1 - parameter value delta1
%       delta2 - parameter value delta2
%   Output:
%       f - function value
% Author: Daphne Nesenberend, Alexey Kazarnikov
% Date: April 2025
function f = f(K,eta,delta1,delta2)


if K <= delta1-1  
    f = eta;
elseif delta1 -1 <= K && K <= delta1 + delta2 -1 
    f = -eta*(K-(delta1-1))^2/(2*delta2*(1-2*delta1-delta2))+eta;
elseif  delta1 + delta2 -1 <= K && K <= - delta1 - delta2
    f = -eta/(1-2*delta1-delta2)*(K+delta1+0.5*delta2);
elseif  -delta1 -delta2 <= K && K < -delta1
    f = eta*(K+delta1)^2/(2*delta2*(1-2*delta1-delta2));
elseif K >= -delta1
    f = 0;
end

end

