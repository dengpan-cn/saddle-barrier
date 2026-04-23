function [newsimInfo, iter] = up_search(simInfo, m, k)
% up_search.m - Perform an upward search from an m-index saddle to a
%               k-index saddle.
%
% Author: Yuchen Xie
% This file is adapted and modified from earlier code by Jianyuan Yin.
% Copyright (c) 2024-2026 Yuchen Xie.
%
% Inputs:
% simInfo:      current state
% m:            initial saddle index
% k:            target saddle index, with k > m
%
% Outputs:
% newsimInfo:   relaxed configuration returned by HiSD
% iter:         iteration counter reported by hisd_sirqit

%% Settings.
options = hisdoptions(1);
options.innpfunc = @(x, y)(x' * y);
options.maxiter = 100000;
options.k =  k ;
simInfo = jampel_calcForce(simInfo);
simInfo = jampel_calcHessian(simInfo);
[blockVectorX, ~] = eigs(simInfo.hessian, k, 'smallestreal');
options.betau = 0.1; options.betat = Inf; options.dt = 0.01;
options.gammamin = 1; options.gammamax = Inf; 
options.outputX = 0; options.outputp = 0; options.outputd = 0;

%% Initialize eigen vectors.
i = randi(k-m);
u = blockVectorX(:, i+m);
V = blockVectorX(:, 1:options.k);
options.V = V;

%% Disturb the saddle to obtain the initial value.
x = reshape(simInfo.xyz', [], 1);
x = x + 10^(-8)*u*rand();
tmp_xyz = reshape(x, 3, []);
simInfo.xyz = tmp_xyz';
simInfo = wrapPos(simInfo);

%% Upward search.
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
