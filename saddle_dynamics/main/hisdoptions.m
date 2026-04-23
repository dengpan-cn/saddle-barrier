function options = hisdoptions(x)
% hisdoptions.m - Build a default option structure for the HiSD solver.
%
% Author: Yuchen Xie
% This file is adapted and modified from earlier code by Jianyuan Yin.
% Copyright (c) 2024-2026 Yuchen Xie.
%
% Input:
% x         relaxation parameter, typically x >= 1
%
% Common fields in the returned options structure:
% V         initial value of V; n-by-k
% innpfunc  inner product function returning one value per column
% k         index of the target saddle point
% dt        initial step size; if dt = 0, only the V iteration is used
% l         initial dimer length
% minl      minimal dimer length
% epsf      stopping tolerance for |F(x)|
% maxiter   maximal iteration count
% betat     maximal beta
% betau     minimal beta
% gammamax  maximal gamma
% gammamin  minimal gamma
% outputX   save iteration data every outputX steps
% outputp   print iteration data every outputp steps
% outputd   draw iteration data every outputd steps
% draws     drawing command executed when outputd is enabled
% tau       maximal displacement magnitude used in step-size control
% adding    additional command executed at each iteration
% 
% Version:	2022.09.08
% Create:	2018.08.07
% Coder:	Jianyuan Yin

options = struct('k', 1, 'dt', 0.1, 'l', 1e-6, 'minl', 1e-6, ...
    'epsf', 3e-13, 'maxiter', 1e5, ...
    'betat', max(x, 1), 'betau', min(1, 1/x), 'gammamax', max(x, 1), 'gammamin', min(1, 1/x), ...
    'outputX', 1, 'outputp', 1, 'outputd', 1, 'draws', '', ...
    'tau', 0.5, 'precd', @(p)p);

if x == Inf
    options.betat = Inf;
    options.betau = 0.01;
    options.gammammax = Inf;
    options.gammammin = 1;
end
end
