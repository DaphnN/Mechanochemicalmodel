%cycle_index implements an index shift for the setting of periodic boundary
%conditions
%   Input:
%       k - index value
%       N - number of discretisation points
%   Output:
%       result_k - corrected index value
% Author: Daphne Nesenberend, Alexey Kazarnikov
% Date: April 2025

function [result_k] = cycle_index(k,N)
if (k > 2*N) || (k<=-N)
    error('Index was outside acceptable region!');
end

if (k < 1) 
    result_k = N + k;
elseif k > N 
    result_k = k - N;
else
    result_k = k;
end

end

