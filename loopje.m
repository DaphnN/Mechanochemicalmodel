% loopje.m is a short programme used to loop over parameter values and run
% Numericaltest2.m
%
% Author: Daphne Nesenberend, Alexey Kazarnikov
% Date: April 2025

clc
clear all
close all


par.D = 0.0001;%Diffusion constant of the morphogen
par.initialN = 256;
par.B = 1;
par.n = 1; %domain scaling
par.L1 = 3 ;
par.L = par.n * par.L1;%domain lenght
par.show = false; %plotting intermediate solutions
par.continuousf = true;
par.ampl =1.5;%parameter for controling length C
par.factorN = 1; %scaling of N
par.N = par.factorN*par.initialN; %amount of discretization points
par.delta1 =0.01;
par.delta2=par.delta1;
par.beta=2;
par.alpha=1;
par.eta =1;

%loop can be adapted to desired parameter
for D = [0.0001]
    par.D=D;
    NumericalTest2(par)
end