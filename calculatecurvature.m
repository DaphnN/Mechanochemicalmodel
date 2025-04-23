% calculatecurvature computes the curvature and its derivatives
%   Input:
%       w - state vector (height, morphogen)
%       par - structure array containing parameters (see loopje.m)
%
%   Output:
%       K - vector with curvature, u, second derivative of the curvature, 
%       first derivative of the curvature and v
%
% Author: Daphne Nesenberend, Alexey Kazarnikov
% Date: April 2025
function [K] = calculatecurvature(w,par)
N = par.N;
beta = par.beta;
h = par.L / par.N;

if (mod(length(w),2)~=0)
    error('Length of input vector w must be an even number!');
end

% array extraction
u = w(1:N);
phi = w(N+1:end);


% tissue curvature
k_fun = @(du,d2u) d2u .* (1 + du.^2).^(-3/2);

% first spatial derivative of tissue curvature
dk_fun = @(du,d2u,d3u) d3u .* (1 + du.^2).^(-3/2) ...
    - 3 * du .* (d2u).^2 .* (1 + du.^2).^(-5/2);


% preferred local curvature
k_hat_fun = @(phi) - beta * phi;

% first spatial derivative of the preferred local curvature
dk_hat_fun = @(dphi) - beta * dphi;

K = zeros(6,N);

for k = 1:N    
        k_p = cycle_index(k+1,N);
        k_pp = cycle_index(k+2,N);
        k_m = cycle_index(k-1,N);
        k_mm = cycle_index(k-2,N);

        du = (u(k) - u(k_m))/h;
        d2u = (u(k_p) - 2*u(k) + u(k_m))/h^2;
        d3u = (u(k_pp) - 2*u(k_p) + 2*u(k_m) - u(k_mm))/(2*h^3);
        dphi = (phi(k) - phi(k_m))/h;

    K(1,k) = k_fun(du,d2u);%K
    K(2,k) = k_fun(du,d2u)-k_hat_fun(phi(k));% u
    K(3,k) = (dk_fun(du,d2u,d3u)-dk_hat_fun(dphi))/sqrt(1+du^2); 
end

for k = 1:N    
    k_p = cycle_index(k+1,N);
    k_m = cycle_index(k-1,N);
    K(4,k) = (K(1,k_p) - K(1,k_m))/(2*h);% first derivative of K
    K(5,k) = (K(1,k_p)-2*K(1,k)+K(1,k_m))/h^2;% second derivative of K 
    K(6,k) = (K(2,k)-K(2,k_m))/h; %v
end

end