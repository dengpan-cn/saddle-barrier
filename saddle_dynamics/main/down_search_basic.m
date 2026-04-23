function [newsimInfo, iter] = down_search_basic(simInfo, m, k, ind, ii)
% down_search_basic.m - Perform a downward search from an m-index saddle
%                       to a k-index saddle along a selected unstable mode.
%
% Author: Yuchen Xie
% This file is adapted and modified from earlier code by Jianyuan Yin.
% Copyright (c) 2024-2026 Yuchen Xie.
%
% Inputs:
% simInfo:      current state
% m:            initial saddle index
% k:            target saddle index, with k < m
% ind:          selected unstable eigendirection for the perturbation
% ii:           perturbation sign, 1 for positive and -1 for negative
%
% Outputs:
% newsimInfo:   relaxed configuration returned by HiSD
% iter:         iteration counter reported by hisd_sirqit

%% Settings.
options = hisdoptions(1);
options.innpfunc = @(x, y)(x' * y);
options.maxiter = 100000;
options.k = k;
options.betau = 0.1; options.betat = Inf; options.dt = 0.01;
options.gammamin = 1; options.gammamax = Inf; 
options.outputX = 0; options.outputp = 0; options.outputd = 0;

%% Calculate the disturbance direction.
simInfo = jampel_calcForce(simInfo);
simInfo = jampel_calcHessian(simInfo);
[blockVectorX, ~] = eigs(simInfo.hessian, m, 'smallestreal');
u = blockVectorX(:, ind);

%% Calculate the initial eigen vectors.
M = 1:m;
if k > 0
    targ = nchoosek(M, k);
    [a, ~] = size(targ);
    V = blockVectorX(:, targ(randi(a), :));
else
    V = blockVectorX(:, 1:options.k);
end
options.V = V;

%% Disturb the saddle to obtain the initial value.
x = reshape(simInfo.xyz', [], 1);
r = rand();
if ii == 1
    x = x + 10^(-8)*u*r;
else
    x = x - 10^(-8)*u*r;
end
tmp_xyz = reshape(x, 3, []);
simInfo.xyz = tmp_xyz';
simInfo = wrapPos(simInfo);

%% Downward search.
[newsimInfo, iter, ~] = hisd_sirqit(simInfo, options);
newsimInfo = wrapPos(newsimInfo);
newsimInfo = jampel_calcForce(newsimInfo);
newsimInfo = jampel_calcHessian(newsimInfo);

end


%% Wrap particle positions using periodic boundary conditions.
function simInfo = wrapPos(simInfo)
simInfo.xyz(:, 1) = simInfo.xyz(:, 1) - round(simInfo.xyz(:, 1)/simInfo.simBoxLen(1)) * simInfo.simBoxLen(1);
simInfo.xyz(:, 2) = simInfo.xyz(:, 2) - round(simInfo.xyz(:, 2)/simInfo.simBoxLen(2)) * simInfo.simBoxLen(2);
simInfo.xyz(:, 3) = simInfo.xyz(:, 3) - round(simInfo.xyz(:, 3)/simInfo.simBoxLen(3)) * simInfo.simBoxLen(3);
end
