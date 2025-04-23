% NumericalTest2 simulates the evolution equations for the 
% 1D mechanochemical model using the implicit TR-BDF2 method. 
% The model is discretised by the Method of Lines,
% thus spatial derivatives are first discretised by a finite difference
% scheme and next the resulting ODE system is solved numerically. Implicit
% method is used to properly control the global arc length integral
% constraint.
%
%   Input:
%     par - structure with model parameter values (see loopje.m)
%     
%   Output:
%     Folder with plots of initial conditions, steady states and
%     variables.mat
% Author: Daphne Nesenberend, Alexey Kazarnikov
% Date: April 2025

function Num2 = NumericalTest2(par)

L = par.L1; % length of the spatial domain
n = par.n; %length change factor (factor that we can also use to change i.c.)
nL = n*L;
N = par.N; % number of discretisation points
h = nL/N; % spatial step, 
xr1 = linspace(0,nL,N+1); % spatial points, where solution will be computed
xr = xr1(1:N);


% defining initial data for simulation
ampl = 0.001 + 0.004 * [par.ampl*0.9646    par.ampl*1    0.8895    0.3416];
freq = 1 + [1     2     1     3];
shift = 0.01 + 0.01 * [0.8    0.3    0.7    0.1];

u0_fun = @(x) 1 ...
    + ampl(1)/par.factorN*sin(freq(1)*2*pi*par.factorN*x/L - shift(1)) ...
    + ampl(2)/par.factorN*sin(freq(2)*2*pi*par.factorN*x/L + shift(2));
phi0_fun = @(x) 0.01 ...
    + ampl(3)/par.factorN*sin(freq(3)*2*pi*par.factorN*x/L + shift(3)) ...
    + ampl(4)/par.factorN*sin(freq(4)*2*pi*par.factorN*x/L - shift(4));

% defining vector of initial conditions
u0 = u0_fun(xr);
phi0 = phi0_fun(xr);
w = [u0'; phi0'; 0]; % the last component stands for Lagrange multiplier
C1=0;
for i = 1:N-1
    C = C1+sqrt(h^2+(u0(i)-u0(i+1))^2);
    C1=C;
end
Clength=C+sqrt(h^2+(u0(1)-u0(N))^2)

% defining model parameters
A=par.eta*par.beta/par.alpha

% create folder for output data
foldername = sprintf('del1=%f, del2=%f, A=%.2f, D=%f, N=%.0f, L=%.1f, C=%f, beta=%.1f, alpha=%.1f',par.delta1, par.delta2, par.beta*par.eta/par.alpha, par.D, par.N, par.L, Clength, par.beta, par.alpha);
mkdir(foldername);

% plotting initial data
figure(1)
clf
hold on
plot(xr,u0,'r.-')
plot(xr,phi0,'b.-')
hold off
xlabel('x')
title('Initial data')
legend('h(x,0)','\phi(x,0)')
drawnow;
saveas(figure(1),[foldername '/Initial_data.fig'])



%Also plot the initial curvature 
K = calculatecurvature(w(1:end-1),par);

figure(2)
    clf
    plot(xr,K(1,:))
    title(sprintf('lambda=%f, A=%.2f, D=%f, N=%.0f, L=%.1f, C=%f, beta=%.1f, eta=%.1f',w(end), par.beta*par.eta/par.alpha, par.D, par.N, par.L, Clength, par.beta, par.eta));
    xlabel('x')
    lgd2 = legend('K');
    set(lgd2,'FontSize',14);
    set(gca, 'FontSize', 14)
    drawnow;
saveas(figure(2),[foldername '/Initial_data_curvature.fig'])


% defining function for global arc length
S = @(u) trapz(xr,sqrt(1+((u - [u(end); u(1:end-1)])./h).^2));

% defining numerical solver parameters
t = 0; % initial time
dt = 1e-5; % initial time step
wprev = w; % previous point of the trajectory
gamma = 2 - sqrt(2); % TR-BDF2 constant
itmax = 50; % maximal allowed number of iterations for Newton's method, used
% to minimise non-linear residual functions inside the numerical method
omega = 1; % relaxation constant
abstol = 1e-5; % absolute tolerance
reltol = 1e-5; % relative tolerance
relerr = 1e-3; % local error
jac_update_interval = 50; % number of iterations during which the Jacobian 
% matrix is kept constant
dFmax = 1e-5; % a threshold for the norm of r.h.s. time derivative used to 
% control the convergence to the steady state pattern and terminate the
% iterations. 
max_step = 5000; % maximal number of time steps
p = 2; % Lp-norm power, used to evaluate vector norms

% this variable is used to collect the statistics of the numerical method.
% Comment this line for better performance.
data_bank = []; 

rhs_fun = @(y) evaluate_rhs_mc_periodic(y(1:end-1),y(end),par);

for k=1:max_step
    % defining method constants, which depend on time step value
    alpha0 = (1-gamma)/(gamma*dt);
    alpha1 = 1/(gamma*(1-gamma)*dt);
    alpha2 = (2-gamma)/((1-gamma)*dt);
    
    % defining residual functions for numerical method (see documentation for
    % TR-BDF2 for details)
    F1 = @(wm,wprev,par) [...
        wprev(1:end-1) + gamma*dt/2 *...
        (...
            rhs_fun(wprev) +...
            rhs_fun(wm)...
        ) - wm(1:end-1);
        globalarclength(wm(1:par.N),par.N,h) - globalarclength(wprev(1:par.N),par.N,h);...
        ];
    G = @(wm,par) [ ...
        gamma*dt/2 * rhs_fun(wm) - wm(1:end-1); ...
        globalarclength(wm(1:par.N),par.N,h) ...
        ];
    
    F2 = @(wm,wint,wprev,par) [...
        alpha1/alpha2*wint(1:end-1) - alpha0/alpha2*wprev(1:end-1) + ...
        (1/alpha2)*rhs_fun(wm) - ...
        wm(1:end-1);...
        globalarclength(wm(1:par.N),par.N,h)-globalarclength(wprev(1:par.N),par.N,h)];
    
    reject_flag = 0; % step rejection flag
    
    % intermediate steps of the method
    wm1 = wprev;
    wm2 = wprev;
    
    finish1 = false; % stage 1 finishing flag
    finish2 = false; % stage 2 finishing flag
    
    % both stages of the numerical method are computed in the same loop by
    % simplified Newton's iterations, in order to gain more performance.
    for it=1:itmax    
        % we evaluate the Jacobian matrix only once per specified number of
        % iterations
        if mod(it,jac_update_interval) == 1
            Jm = evaluate_jacobian_num(@(w) G(w,par),wm2,1e-8);
            
            if exist('data_bank','var')
                % eigenvalues of the Jacobian matrix
                Jm_eig = eig(Jm);
                eig_r = real(Jm_eig);
                eig_im = imag(Jm_eig);
                
                % closest real part to zero
                [~,min_idx] = min(abs(eig_r));
                
                % matrix condition number
                Jm_cond = cond(Jm);
            end
        end
        
        if ~finish1 % first stage
            Fwm1 = F1(wm1,wprev,par);
            dw1 = Jm \ -Fwm1;
            wnew1 = wm1 + omega * dw1;
            
            % divergence control
            Fwnew1 = F1(wnew1,wprev,par);
            if vecnorm(Fwm1,p) <= vecnorm(Fwnew1,p)
                reject_flag = true;
                break;
            end
            
            % convergence control
            if it > 1 && ...
                vecnorm(Fwnew1,p) < abstol && ...
                vecnorm(wm1-wnew1,p) < reltol
            
                finish1 = true;
            end
            
            wm1 = wnew1;
        end
        
        % second stage
        Fwm2 = F2(wm2,wm1,wprev,par);                            
        dw2 = Jm \ -Fwm2;
        wnew2 = wm2 + omega * dw2;    
        Fwnew2 = F2(wnew2,wnew1,wprev,par);

        % convergence control
        if it > 1 && finish1 && vecnorm(Fwnew2,p) < abstol && vecnorm(wm2-wnew2,p) < reltol
            break;
        end
        
        % Note that due to joint computing of both stages of the method,
        % second part of the iterations can diverge at the beginning, thus
        % divergence is not tested here.

        wm2 = wnew2;

        % can be used for debugging purposes
        %disp(norm(Fwm1))
        %disp(norm(Fwm2))   
    end
 
    % collecting statistics for intermediate step
    if exist('data_bank','var')
        didx = 1;
        col_num = size(data_bank,2)+1;
        data_bank(didx,col_num) = Jm_cond; didx = didx+1;
        data_bank(didx,col_num) = norm(Fwnew1); didx = didx+1;
        data_bank(didx,col_num) = reject_flag || it == itmax; didx = didx+1;
        data_bank(didx,col_num) = norm(wm1-wnew1); didx = didx+1;
        data_bank(didx,col_num) = eig_r(min_idx); didx = didx+1;
        data_bank(didx,col_num) = eig_im(min_idx); didx = didx+1;
        data_bank(didx,col_num) = it; didx = didx+1; 
        data_bank(didx,col_num) = norm(Fwnew2); didx = didx+1;
        data_bank(didx,col_num) = norm(wm2-wnew2); didx = didx+1;
        data_bank(didx,col_num) = t; didx = didx+1;
        data_bank(didx,col_num) = dt; didx = didx+1;
    end
    
    % if divergence was detected or maximal number of iterations was
    % reached, current step will be rejected
    if reject_flag || it == itmax
        dt = 0.5*dt;
        continue;
    end
    
    wint = wnew1;
    wnew = wnew2;
    
    % time step size control (see documentation for TR-BDF2 for details)  
    rhs = rhs_fun(wprev);  
    tau = (3*gamma^2 - 4*gamma + 2)/(6*(1-gamma)^2) * ...
        (wnew(1:end-1) ...
        - 1/gamma^2*wint(1:end-1) ...
        + (1-gamma^2)/(gamma^2)*wprev(1:end-1) + ...
        dt*(1-gamma)/gamma*rhs ...
        ); 
    taur = abs(tau ./ wnew(1:end-1));
    t_opt = dt*(relerr/norm(taur))^(1/3);
    t_opt = min(t_opt, 2 * dt); % prevent too rapid time step growth
    
    % advancing to the next time step
    t = t + dt;
    wprev = wnew; 
    dt = 0.9 * t_opt;
    
    % plotting intermediate solution
    if par.show
         K = calculatecurvature(wprev(1:end-1),par);
        figure(3)
        clf
        plot(xr,K(1,:))
        title(sprintf('lambda=%f, A=%.2f, D=%f, N=%.0f, L=%.1f, C=%f, beta=%.1f, eta=%.1f',wprev(end), par.beta*par.eta/par.alpha, par.D, par.N, par.L, Clength, par.beta, par.eta));
        xlabel('x')
        lgd2 = legend('K');
        set(lgd2,'FontSize',14);
        set(gca, 'FontSize', 14)
        drawnow;

        figure(4)
        clf
        hold on
        plot(xr,wprev(1:N),'r.-')
        plot(xr,wprev(N+1:end-1),'b.-')
        hold off
        title(sprintf('t=%f, S=%f, dt=%f',t, S(wprev(1:N)), dt))   
        xlabel('x')
        lgd = legend('h(x,t)','\phi(x,t)');
        set(lgd,'FontSize',14);
        set(gca, 'FontSize', 14)
        drawnow;
    end

    % convergence control
    time_deriv = rhs_fun(wprev);
    norm_deriv = vecnorm(time_deriv(N+1:end),p);
    if norm_deriv < dFmax
        break
    end   
    Wdata(k,:)=wprev;
end
outputdhdt=evaluate_rhs_mc_periodic(wprev(1:end-1),wprev(end),par);
%produce final figures:

K = calculatecurvature(wprev(1:end-1),par);

figure(3)
    clf
    plot(xr,K(1,:))
    title(sprintf('lambda=%f, A=%.2f, D=%f, N=%.0f, L=%.1f, C=%f, beta=%.1f, eta=%.1f',wprev(end), par.beta*par.eta/par.alpha, par.D, par.N, par.L, Clength, par.beta, par.eta));
    xlabel('x')
    lgd2 = legend('K');
    set(lgd2,'FontSize',14);
    set(gca, 'FontSize', 14)
    drawnow;
saveas(figure(3),[foldername '/Final_curvature.fig']) 

    
 figure(4)
    clf
    hold on
    plot(xr,wprev(1:N),'r.-')
    plot(xr,wprev(N+1:end-1),'b.-')
    hold off
    title(sprintf('dhdt=%f, lambda=%f, A=%.2f, D=%f, N=%.0f, L=%.1f, C=%f, beta=%.1f, eta=%.1f',outputdhdt(1), wprev(end), par.beta*par.eta/par.alpha, par.D, par.N, par.L, Clength, par.beta, par.eta));
    xlabel('x')
    lgd = legend('h(x,t)','\phi(x,t)');
    set(lgd,'FontSize',14);
saveas(figure(4),[foldername '/Final_h_phi']);
    
 figure(5)
    clf
    plot3(par.beta*wprev(N+1:end-1),K(6,:),K(2,:),'-o') 
    xlabel('\psi')
    ylabel('v')
    zlabel('u')
    title(sprintf('dhdt=%f, lambda=%f, A=%.2f, D=%f, N=%.0f, L=%.1f, C=%f, beta=%.1f, eta=%.1f',outputdhdt(1), wprev(end), par.beta*par.eta/par.alpha, par.D, par.N, par.L, Clength, par.beta, par.eta));
saveas(figure(5),[foldername '/final_u_v_psi_plane']);


 %save data (in correct folder):
 save([foldername '/variables.mat'])
 
disp('All done!')

