function flag = jampel_compareConf(s1, s2)
% jampel_compareConf.m - Check whether two configurations are identical.
%
% Author: Yuchen Xie
% Copyright (c) 2024-2026 Yuchen Xie.
%
% Inputs:
% s1, s2:        the two configurations
%
% Output:
% flag:          true if the configurations are equivalent, false otherwise

%% Settings.
n1 = size(s1.xyz, 1);
if log2(n1) ~= floor(log2(n1))
    s1.xyz = [s1.xyz; [0, 0, 0]];
end
n2 = size(s2.xyz, 1);
if log2(n2) ~= floor(log2(n2))
    s2.xyz = [s2.xyz; [0, 0, 0]];
end

if any(size(s1.ratIdx) ~= size(s2.ratIdx))
    flag = false;
    return;
end

% s1 and s2 has the same number of rattlers
nr = size(s1.ratIdx, 2);
aPara = 0.3 * (2.0 * s1.radius);
[N, ~] = size(s1.xyz);
qcut = 1e-2;

%% Calculate dist matrix.
distMat = zeros([N, N]);
for iatom = 1:N
    dr = [];
    dr(:, 1) = s1.xyz(iatom, 1) - s2.xyz(:, 1);
    dr(:, 2) = s1.xyz(iatom, 2) - s2.xyz(:, 2);
    dr(:, 3) = s1.xyz(iatom, 3) - s2.xyz(:, 3);
    % periodic boundary condition
    dr(:, 1) = dr(:, 1) - round(dr(:, 1)/s1.simBoxLen(1)) * s1.simBoxLen(1);
    dr(:, 2) = dr(:, 2) - round(dr(:, 2)/s1.simBoxLen(2)) * s1.simBoxLen(2);
    dr(:, 3) = dr(:, 3) - round(dr(:, 3)/s1.simBoxLen(3)) * s1.simBoxLen(3);
    
    distMat(iatom, :) = (dr(:, 1).^2 + dr(:, 2).^2 + dr(:, 3).^2).^0.5;
end
distMat = exp(-distMat/aPara);
% scale weights of rattlers
distMat(s1.ratIdx, s2.ratIdx) = distMat(s1.ratIdx, s2.ratIdx) / (0.367879 * 0.367879);
distMat(s1.ratIdx, :) = distMat(s1.ratIdx, :) * 0.367879;
distMat(:, s2.ratIdx) = distMat(:, s2.ratIdx) * 0.367879;
distMat(distMat<qcut) = 0;

%% Max matching.
[m, ur, uc] = matchpairs(distMat, -N, 'max');     % cheaper
% any un-matching means difference
if ~(isempty(ur) && isempty(uc))
    flag = false;
    return;
end

%% Checking.
% rattler to rattler?
for ith = 1:nr
    tmp = m(m(:, 1) == s1.ratIdx(ith), :);
    if ~any(s2.ratIdx == tmp(2))
        flag = false;
        return;
    end
end

lnr = 1:N;
lnr(s1.ratIdx) = 0;

dr = s1.xyz(m(lnr > 0, 1), :) - s2.xyz(m(lnr > 0, 2), :);
dr(:, 3) = dr(:, 3) - round(dr(:, 3)/s1.simBoxLen(3)) * s1.simBoxLen(3);
dr(:, 2) = dr(:, 2) - round(dr(:, 2)/s1.simBoxLen(2)) * s1.simBoxLen(2);
dr(:, 1) = dr(:, 1) - round(dr(:, 1)/s1.simBoxLen(1)) * s1.simBoxLen(1);
% coinciding the center of these two configurations
dr(:, 1) = dr(:, 1) - mean(dr(:, 1));
dr(:, 2) = dr(:, 2) - mean(dr(:, 2));
dr(:, 3) = dr(:, 3) - mean(dr(:, 3));

dr = sqrt(dr(:, 1).^2 + dr(:, 2).^2 + dr(:, 3).^2);
if all(dr <= 1E-9 * 2.0 * s1.radius)
    flag = true;
else
    flag = false;
end

end


