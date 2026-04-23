function [newsimInfo, iter, output] = hisd_sirqit(simInfo, options)
% hisd_sirqit.m - Run the HiSD solver with SIRQIT-style direction updates.
%
% Author: Yuchen Xie
% This file is adapted and modified from earlier code by Jianyuan Yin.
% Copyright (c) 2024-2026 Yuchen Xie.
%
% Input:
% simInfo:              current state
%
% options:
% innpfunc:             inner product function handle with two inputs (recommended)
% k:                    index of hisd (recommended)
% V:                    initial direction vectors, should be k column vectors (recommended)
% dt:                   first stepsize (optional, 0.01)
% l:                    dimer length (optional, 1e-6)
% minl:                 minimal dimer length (optional, 1e-6)
% epsf:                 tolerance of gradient (optional, 1e-8)
% maxiter:              maximal iteration steps (optional, 1e4)
% betat, betau:         max/minimal beta (optional, 1, 1)
% gammamax, gammamin:   max/minimal gamma (optional, 1, 1)
% outputX:              save iteration messages every outputX steps (optional, 0 not save)
% outputp:              print iteration message every outputp steps (optional, 0 not print)
% outputd, draws:       draw iteration message using draws every outputd steps (optional, 0 not draw)
% tau:                  max distance for every iteration (optional, 0.1)
% adding:               extra evaluation for each iteration (optional, '')
%
% Outputs:
% newsimInfo:           updated configuration after the solver stops
% iter:                 solver status / iteration count used by this codebase
% output:               diagnostic data collected during the run
%
% Version:	2022.09.08
% Create:	2018.08.07
% Coder:	Jianyuan Yin

%% Reading parameters.
x = reshape(simInfo.xyz', [], 1);
n = length(x);          % n is the number of variables

if isfield(options, 'innpfunc')
    innp = @(x, y)options.innpfunc(x, y);
    nor = @(x)sqrt(innp(x, x));
else
    disp('HiSD: using the standard inner product.');
    innp = @(x, y)(x' * y);
    nor = @(x)sqrt(x' * x);
end

if isfield(options, 'k')
    if isfield(options, 'V')
        if options.k ~= size(options.V, 2)
            options.k = size(options.V, 2);
            disp(['HiSD: k is reset as ' num2str(options.k) ', the column - number of V.']);
        end
    end
elseif isfield(options, 'V')
    options.k = size(options.V, 2);
else
    options.k = 1;
    disp('HiSD: k is reset as 1.');
end
k = options.k;

if isfield(options, 'V')
    if rank(options.V) < k
        disp('HiSD: V is random vectors.');
        V = randn(n, k);
    else
        V = options.V;
    end
else
    disp('HiSD: V is random vectors.');
    V = randn(n, k);
end

% Gram-Schmidt orthogonalization
for i = 1 : k
    V(:, i) = V(:, i) - V(:, 1:i-1) * innp(V(:, 1:i-1), V(:, i));
    V(:, i) = V(:, i)/nor(V(:, i));
end

if isfield(options, 'dt')
    dt = options.dt;
else
    disp('HiSD: dt = 0.01.');
    dt = 0.01;
end

if isfield(options, 'l')
    l = options.l;
else
    l = 1e-6;
end

if isfield(options, 'minl')
    minl = options.minl;
else
    minl = 1e-6;
end

if isfield(options, 'epsf')
    epsf = options.epsf;
else
    disp('HiSD: epsf = 1e-8.');
    epsf = 1e-8;
end

if isfield(options, 'maxiter')
    maxiter = options.maxiter;
else
    maxiter = 1e4;
end

if isfield(options, 'betat')
    betat = options.betat;
else
    betat = 1;
end
if isfield(options, 'betau')
    betau = options.betau;
else
    betau = 1;
end

if ~isfield(options, 'outputX')
    options.outputX = 0;
end

if options.outputX ~= 0
    X = zeros(k+2+n, ceil(maxiter/options.outputX));
    storing = 'X(:,ceil(iter/options.outputX))=[nor(f);alpha;beta;x];';
else
    storing = '';
end

if ~isfield(options, 'outputp')
    options.outputp = 0;
end

if ~isfield(options, 'outputd')
    options.outputd = 0;
end

if options.outputd ~= 0 && isfield(options, 'draws')
    drawing = options.draws;
else
    drawing = '';
end

if isfield(options, 'tau')
    tau = options.tau;
else
    tau = 0.1;
end

if isfield(options, 'adding')
    adding = options.adding;
else
    adding = '';
end

if isfield(options, 'gammamax')
    gammamax = options.gammamax;
else
    gammamax = 1;
end
if isfield(options, 'gammamin')
    gammamin = options.gammamin;
else
    gammamin = 1;
end

%% Setting environment.
alpha (1 : k, 1) = 0;
beta = dt;
gamma (1 : k) = dt;
dV = zeros(n, k);
F = @ (x) Force(x, simInfo);
H = @ (xx, vv, ll) - (F(xx + ll * vv) - F(xx - ll * vv)) / (2 * ll);
f = F(x);

flag = 0;

%% First iteration.
% iter = 1;

% update x
g = f - 2 * V * innp(V, f);
gp = g;
Dx = g * beta;
oldx = x;
x = x + Dx;
momentum = x - oldx;

% periodic boundary condition
momentum_matrix = reshape(momentum, 3, []);
momentum_matrix = momentum_matrix';
momentum_matrix(:, 1) = momentum_matrix(:, 1) - round(momentum_matrix(:, 1)/simInfo.simBoxLen(1)) * simInfo.simBoxLen(1);
momentum_matrix(:, 2) = momentum_matrix(:, 2) - round(momentum_matrix(:, 2)/simInfo.simBoxLen(2)) * simInfo.simBoxLen(2);
momentum_matrix(:, 3) = momentum_matrix(:, 3) - round(momentum_matrix(:, 3)/simInfo.simBoxLen(3)) * simInfo.simBoxLen(3);
momentum = reshape(momentum_matrix', [], 1);

for i = 1 : k
    ui = H(x, V(:, i), l);
    alpha(i) = innp(V(:, i), ui);
    dV(:, i)= -ui + alpha(i) * V(:, i) + V(:, 1:i-1) * (2 * innp(V(:, 1:i-1), ui));
end
dVp = dV;

% updata V
Vp = V;
V = V + dV .* gamma;
for i = 1 : k
    V(:, i) = V(:, i) - V(:, 1:i-1) * innp(V(:, 1:i-1), V(:, i));
    V(:, i) = V(:, i)/nor(V(:, i));
end
DV = V - Vp;
% shrink dimer
l = max(l/(1+dt), minl);

f = F(x);
norf = nor(f);
eval(storing);

%% Following BB iterations.
iter = 2;
while iter <= maxiter
    g = f - 2 * V * innp(V, f);
    Dg = g - gp;
    gp = g;
    
    beta = abs (innp(Dx, Dg) / innp(Dg, Dg));
    beta = min (beta, betat * dt);
    beta = max (beta, betau * dt);
    if norf*beta > tau
        beta = tau/norf;
    end
    
    % heavy ball parameters
    if norf < 1e-10 && ~flag
        simInfo = jampel_calcHessian(simInfo);
        m = abs(eigs(simInfo.hessian, 1, 1e-8));
        M = eigs(simInfo.hessian, 1, 'largestabs');
        s = jampel_cvtRun2ArXiv(simInfo);
        sl = length(s.ratIdx);
        HBgamma = (1 - 2 / (sqrt(M / m) + 1)) ^ 2;
        HBbeta = (2 / (sqrt(M) + sqrt(m))) ^ 2;
        flag = 1;
    end
    
    % rattlers
    if flag && sl > 0
        flag = 0;
        fprintf('m = %i \n', m);
        break;
    end
    
    % heavy ball
    if norf < 1e-11
        gamma0 = HBgamma;
        beta = HBbeta;
    else
        gamma0 = 0;
    end
    
    % update x
    Dx = g * beta;
    oldx = x;
    x = x + Dx + gamma0 * momentum;
    
    momentum = x - oldx;
    % periodic boundary condition
    momentum_matrix = reshape(momentum, 3, []);
    momentum_matrix = momentum_matrix';
    momentum_matrix(:, 1) = momentum_matrix(:, 1) - round(momentum_matrix(:, 1)/simInfo.simBoxLen(1)) * simInfo.simBoxLen(1);
    momentum_matrix(:, 2) = momentum_matrix(:, 2) - round(momentum_matrix(:, 2)/simInfo.simBoxLen(2)) * simInfo.simBoxLen(2);
    momentum_matrix(:, 3) = momentum_matrix(:, 3) - round(momentum_matrix(:, 3)/simInfo.simBoxLen(3)) * simInfo.simBoxLen(3);
    momentum = reshape(momentum_matrix', [], 1);
    
    tmp_xyz = reshape(x, 3, []);
    simInfo.xyz = tmp_xyz';
    
    for i = 1:k
        ui = H(x, V(:, i), l);
        alpha(i) = innp(V(:, i), ui);
        dV(:, i) = - ui + alpha(i) * V(:, i) + V(:, 1:i-1) * (2 * innp(V(:, 1:i-1), ui));
    end
    Dd = dVp - dV;
    dVp = dV;
    
    for i = 1:k
        gamma(i) = abs(innp(DV(:, i), Dd(:, i)) / innp(Dd(:, i), Dd(:, i)));
        if ~isfinite(gamma(i))
            gamma(i) = 1;
        end
    end
    gamma = min(gamma, gammamax*dt);
    gamma = max(gamma, gammamin*dt);

    % update V
    Vp = V;
    V = V + dV .* gamma;
    for i = 1 : k
        V(:, i) = V(:, i) - V(:, 1:i-1) * innp(V(:, 1:i-1), V(:, i));
        V(:, i) = V(:, i)/nor(V(:, i));
    end
    DV = V - Vp;
    % shrink dimer
    l = max(l/(1+dt), minl);
    
    f = F(x);
    norf = nor(f);
    
    eval(storing);
    if mod(iter, options.outputp) == 0
        disp([num2str(iter, '%d') ': nf=' num2str(norf, '%.4e') '; a=' num2str(alpha', '%.4f ')]);
    end
    if mod(iter, options.outputd) == 0
        eval(drawing);drawnow;
    end
    if norf < epsf && dt > 0
        if options.outputX ~= 0
            X = X(:, 1:ceil(iter/options.outputX));
        end
        break;
    end
    iter = iter + 1;
    eval(adding);
    
    if mod(iter, 1000) == 0
        fprintf('iter = %i, nf = %i, gap = %i \n', iter, norf, norm(momentum, 2));
    end
end

tmp_xyz = reshape(x, 3, []);
simInfo.xyz = tmp_xyz';
simInfo = wrapPos(simInfo);
newsimInfo = simInfo;

if iter == maxiter + 1
    iter = maxiter;
end
if ~flag
    iter = -1;
end

if dt == 0
    x = V;
    output = alpha;
    return
end
output = struct('x', x, 'V', V, 'it', iter, 'al', alpha);
if options.outputX ~= 0
    output.X = X;
end
disp(iter);

end

%% Evaluate forces for a flattened coordinate vector.
function f = Force(x, simInfo)
tmp_xyz = reshape(x, 3, []);
simInfo.xyz = tmp_xyz';
simInfo = jampel_calcForce(simInfo);
f = reshape(simInfo.force', [], 1);
end

%% Wrap particle positions using periodic boundary conditions.
function simInfo = wrapPos(simInfo)
simInfo.xyz(:, 1) = simInfo.xyz(:, 1) - round(simInfo.xyz(:, 1)/simInfo.simBoxLen(1)) * simInfo.simBoxLen(1);
simInfo.xyz(:, 2) = simInfo.xyz(:, 2) - round(simInfo.xyz(:, 2)/simInfo.simBoxLen(2)) * simInfo.simBoxLen(2);
simInfo.xyz(:, 3) = simInfo.xyz(:, 3) - round(simInfo.xyz(:, 3)/simInfo.simBoxLen(3)) * simInfo.simBoxLen(3);
end
