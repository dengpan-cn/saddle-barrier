function dist = jampel_calcDistance(s1, s2)
% jampel_calcDistance.m - Calculate the distance between two configurations.
%
% Author: Yuchen Xie
% Copyright (c) 2024-2026 Yuchen Xie.
%
% Inputs:
% s1, s2:        the two configurations
%
% Output:
% dist:          contact-based configuration distance normalized by diameter

%% Calculate the contact matrix.
cm1 = compContactMatrix(s1);
cm2 = compContactMatrix(s2);

%% Calculate the distance.
lidx = find(cm1 + cm2);
[ii, jj] = ind2sub(size(cm1), lidx);

dr = s1.xyz(ii, :) - s1.xyz(jj, :);
dr(:, 1) = dr(:, 1) - round(dr(:, 1)/s1.simBoxLen(1)) * s1.simBoxLen(1);
dr(:, 2) = dr(:, 2) - round(dr(:, 2)/s1.simBoxLen(2)) * s1.simBoxLen(2);
dr(:, 3) = dr(:, 3) - round(dr(:, 3)/s1.simBoxLen(3)) * s1.simBoxLen(3);
dr1 = dr;

dr = s2.xyz(ii, :) - s2.xyz(jj, :);
dr(:, 1) = dr(:, 1) - round(dr(:, 1)/s2.simBoxLen(1)) * s2.simBoxLen(1);
dr(:, 2) = dr(:, 2) - round(dr(:, 2)/s2.simBoxLen(2)) * s2.simBoxLen(2);
dr(:, 3) = dr(:, 3) - round(dr(:, 3)/s2.simBoxLen(3)) * s2.simBoxLen(3);
dr2 = dr;

dr = diag(cm1(lidx)) * dr1 - diag(cm2(lidx)) * dr2;

dist = sqrt(sum(sum(dr.^2))) / (2.0 * s1.radius);

end


%% Calculate the contact matrix for a configuration.
function contactMat = compContactMatrix(simInfo)
rcutP2 = power(simInfo.radius * 2.0, 2);

n = size(simInfo.xyz,1);
if log2(n) ~= floor(log2(n))
    simInfo.xyz = [simInfo.xyz; [0,0,0]];
end
[N, dim] = size(simInfo.xyz);

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

% remove rattlers
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

contactMat = bondMat;

end
