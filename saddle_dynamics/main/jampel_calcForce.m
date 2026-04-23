function simInfo = jampel_calcForce(simInfo)
% jampel_calcForce.m - Compute forces and related observables for a packing.
%
% Author: Yuchen Xie
% Copyright (c) 2024-2026 Yuchen Xie.
%
% Harmonic / Hertz potential energy
%
% Input / output:
% simInfo:      current state; this function appends force, energy, contact,
%               and pressure-related fields to the structure

%% Settings.
n = size(simInfo.xyz, 1);
if log2(n) ~= floor(log2(n))
    simInfo.xyz = [simInfo.xyz; [0, 0, 0]];
end

rcut = simInfo.radius * 2.0;
N = size(simInfo.xyz, 1);
dim = size(simInfo.xyz, 2);

simInfo.force = zeros([N, dim]);
simInfo.ePair = 0;
simInfo.nContact = 0;
simInfo.pVirTens = zeros([dim, dim]);

%% Calculate force and other properties.
dr = zeros([N-1, dim+1]);
for iatom = 1:N
    dr(:, 1) = simInfo.xyz(iatom, 1) - simInfo.xyz([1:iatom-1, iatom+1:N], 1);
    dr(:, 2) = simInfo.xyz(iatom, 2) - simInfo.xyz([1:iatom-1, iatom+1:N], 2);
    dr(:, 3) = simInfo.xyz(iatom, 3) - simInfo.xyz([1:iatom-1, iatom+1:N], 3);
    % periodic boundary condition
    dr(:, 1) = dr(:, 1) - round(dr(:, 1)/simInfo.simBoxLen(1)) * simInfo.simBoxLen(1);
    dr(:, 2) = dr(:, 2) - round(dr(:, 2)/simInfo.simBoxLen(2)) * simInfo.simBoxLen(2);
    dr(:, 3) = dr(:, 3) - round(dr(:, 3)/simInfo.simBoxLen(3)) * simInfo.simBoxLen(3);
    dr(:, 4) = sqrt(dr(:, 1).^2 + dr(:, 2).^2 + dr(:, 3).^2);

    % for the iatom-th particle, cinfo stores the particles which are contacted with this particle
    cinfo = dr(dr(:, 4) <= rcut, :);
    if isempty(cinfo)
        continue;
    end
    simInfo.nContact = simInfo.nContact + size(cinfo, 1);

    rdivsig = cinfo(:, 4) / rcut;                           % d_{ij}/d

    % Harmonic
    fpair = (1.0 - rdivsig) ./ cinfo(:, 4) / rcut;          % (1-d_{ij}/d)/(d*d_{ij})
    simInfo.ePair = simInfo.ePair + sum(0.5 * (1.0 - rdivsig) .* (1.0 - rdivsig)) / N;

    % Hertz
    % fpair = (1.0 - rdivsig) .^ (3/2) ./ cinfo(:, 4) / rcut;        % (1-d_{ij}/d)^(3/2) / (d*d_{ij})
    % simInfo.ePair = simInfo.ePair + sum(0.4 * (1.0 - rdivsig) .^ (5 / 2)) / N;
    
    simInfo.force(iatom, 1) = sum(fpair .* cinfo(:, 1));
    simInfo.force(iatom, 2) = sum(fpair .* cinfo(:, 2));
    simInfo.force(iatom, 3) = sum(fpair .* cinfo(:, 3));

    tpvir = zeros([dim, dim]);
    for jth = 1:size(cinfo, 1)
        tpvir = tpvir + fpair(jth) * cinfo(jth, 1:3)' * cinfo(jth, 1:3);
    end
    simInfo.pVirTens = simInfo.pVirTens + tpvir;
end

simInfo.ePair = simInfo.ePair / 2.0;                % energy of this configuration
simInfo.nContact = simInfo.nContact / 2;            % number of contacting pairs
simInfo.z = simInfo.nContact * 2.0 / N;
simInfo.pressure = trace(simInfo.pVirTens) / dim;     % pressure of this configuration

%% Fix one particle.
simInfo.force = simInfo.force(1:N-1, :);             % negative gradient
simInfo.xyz = simInfo.xyz(1:N-1, :);

end
