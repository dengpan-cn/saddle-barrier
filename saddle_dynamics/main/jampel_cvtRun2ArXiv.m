function arch = jampel_cvtRun2ArXiv(simInfo)
% jampel_cvtRun2ArXiv.m - Convert a runtime state to a compact archive form.
%
% Author: Yuchen Xie
% Copyright (c) 2024-2026 Yuchen Xie.
%
% Input:
% simInfo:      current state
%
% Output:
% arch:         archived structure containing geometry and topology fields

%% Settings.
arch.xyz = simInfo.xyz;
arch.simBoxLen = simInfo.simBoxLen;
arch.radius = simInfo.radius;
arch.childrenList = [];

% wrap PBC
arch.xyz(:, 1) = arch.xyz(:, 1) - round(arch.xyz(:, 1)/arch.simBoxLen(1)) * arch.simBoxLen(1);
arch.xyz(:, 2) = arch.xyz(:, 2) - round(arch.xyz(:, 2)/arch.simBoxLen(2)) * arch.simBoxLen(2);
arch.xyz(:, 3) = arch.xyz(:, 3) - round(arch.xyz(:, 3)/arch.simBoxLen(3)) * arch.simBoxLen(3);

%% Calculate rattler/buckler.
[arch.ratIdx, arch.buckIdx, ~] = compContactMatrix(arch);

end

%% Compute contact, rattler, and buckler information.
function [ratIdx, buckIdx, contactMat] = compContactMatrix(simInfo)
%% Settings.
rcutP2 = power(simInfo.radius * 2.0, 2);

n = size(simInfo.xyz, 1);
if log2(n) ~= floor(log2(n))
    simInfo.xyz = [simInfo.xyz; [0, 0, 0]];
end
[N, dim] = size(simInfo.xyz);

%% Calculate bond matrix.
bondMat = zeros([N, N]);
for iatom = 1:N-1
    dr = [];
    dr(:, 1) = simInfo.xyz(iatom, 1) - simInfo.xyz(iatom+1:N, 1);
    dr(:, 2) = simInfo.xyz(iatom, 2) - simInfo.xyz(iatom+1:N, 2);
    dr(:, 3) = simInfo.xyz(iatom, 3) - simInfo.xyz(iatom+1:N, 3);
    % periodic boundary condition
    dr(:, 1) = dr(:, 1) - round(dr(:, 1)/simInfo.simBoxLen(1)) * simInfo.simBoxLen(1);
    dr(:, 2) = dr(:, 2) - round(dr(:, 2)/simInfo.simBoxLen(2)) * simInfo.simBoxLen(2);
    dr(:, 3) = dr(:, 3) - round(dr(:, 3)/simInfo.simBoxLen(3)) * simInfo.simBoxLen(3);

    bondMat(iatom, iatom+1:N) = (dr(:, 1).^2 + dr(:, 2).^2 + dr(:, 3).^2 <= rcutP2);
end
bondMat = bondMat + bondMat';

%% Rattlers.
ratIdx = find(sum(bondMat) < dim + 1);
n1 = size(ratIdx, 2);
while (true)
    bondMat(ratIdx, :) = 0;
    bondMat(:, ratIdx) = 0;

    ratIdx = find(sum(bondMat) < dim + 1);
    if n1 == size(ratIdx, 2)
        break;
    end
    n1 = size(ratIdx, 2);
end

if n1 == 0
    ratIdx = [];
end

contactMat = bondMat;

%% Bucklers.
buckIdx = find(sum(bondMat) == dim + 1);
if size(buckIdx, 2) == 0
    buckIdx = [];
end

%% Fix one particle.
simInfo.xyz = simInfo.xyz(1:N-1, :);

end
