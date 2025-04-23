%evaluate_rhs_mc_periodic evaluates the r.h.s. of the mechanochemical model
%in the assumption of periodic boundary conditions
%   Input:
%       w - state vector (height, morphogen)
%       lambda - Lagrange multiplier
%       par - structure with model parameter values (see loopje.m)
%   Output:
%       dwdt - r.h.s. vector
% Author: Daphne Nesenberend, Alexey Kazarnikov
% Date: April 2025

function dwdt = evaluate_rhs_mc_periodic(w,lambda,par)

if (mod(length(w),2)~=0)
    error('Length of input vector w must be an even number!');
end

% structure fields unpacking
N = par.N;
alpha = par.alpha;
beta = par.beta;
D = par.D;
eta = par.eta;
L = par.L;
h = par.L / par.N;

% array extraction
u = w(1:N);
phi = w(N+1:end);

% defining helper functions
% tissue curvature
k_fun = @(du,d2u) d2u .* (1 + du.^2).^(-3/2);

% first spatial derivative of tissue curvature
dk_fun = @(du,d2u,d3u) d3u .* (1 + du.^2).^(-3/2) ...
    - 3 * du .* (d2u).^2 .* (1 + du.^2).^(-5/2);

% second spatial derivative of tissue curvature
d2k_fun = @(du,d2u,d3u,d4u) d4u .* (1 + du.^2).^(-3/2) ...
    - 3 * du .* d2u .* d3u .* (1 + du.^2).^(-5/2) ...
    - 3 * (d2u.^3 + 2 * du .* d2u .* d3u) .* (1 + du.^2).^(-5/2) ...
    + 15 * du.^2 .* d2u.^3 .* (1 + du.^2).^(-7/2);

% preferred local curvature
k_hat_fun = @(phi) - beta * phi;

% first spatial derivative of the preferred local curvature
dk_hat_fun = @(dphi) - beta * dphi;

% second spatial derivative of the preferred local curvature
d2k_hat_fun = @(d2phi) - beta * d2phi;


if par.continuousf
    f_fun = @(du,d2u) f(k_fun(du,d2u),eta,par.delta1,par.delta2); 
else
    f_fun = @(du,d2u)  -eta * max(min(k_fun(du,d2u),0),-1); 
end



dudt = zeros(1,N);
dphidt = zeros(1,N);

for k = 1:N
    k_p = cycle_index(k+1,N);
    k_pp = cycle_index(k+2,N);
    k_m = cycle_index(k-1,N);
    k_mm = cycle_index(k-2,N);
    
    du = (u(k) - u(k_m))/h;
    d2u = (u(k_p) - 2*u(k) + u(k_m))/h^2;
    d3u = (u(k_pp) - 2*u(k_p) + 2*u(k_m) - u(k_mm))/(2*h^3);
    d4u = (u(k_pp) - 4*u(k_p) + 6*u(k) - 4*u(k_m) + u(k_mm))/h^4;
     
    dphi = (phi(k) - phi(k_m))/h;
    d2phi = (phi(k_p) - 2*phi(k) + phi(k_m))/h^2;
    

    dudt(k) = -par.B * ... 
        (...  
(d2k_fun(du,d2u,d3u,d4u) - d2k_hat_fun(d2phi)) ./ (1 + du.^2) ...
        - (dk_fun(du,d2u,d3u) - dk_hat_fun(dphi)).*du.*d2u ./ (1+du.^2).^2 ...
        +(k_fun(du,d2u)-k_hat_fun(phi(k))).*k_fun(du,d2u).^2 ... 
        -0.5.*k_fun(du,d2u).*(k_fun(du,d2u)-k_hat_fun(phi(k))).^2  - lambda * k_fun(du,d2u)); 
    dphidt(k) = D * ( d2phi ./ (1 + du.^2) ...
        - dphi .* du .* d2u ./ (1 + du.^2).^2)...
       - alpha * phi(k) +f_fun(du,d2u);
    
end

dwdt = zeros(1,2*N)';
dwdt(1:N) = dudt;
dwdt(N+1:end) = dphidt;

end

