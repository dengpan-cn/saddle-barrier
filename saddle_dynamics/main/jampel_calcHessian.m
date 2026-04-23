function simInfo = jampel_calcHessian(simInfo)
% jampel_calcHessian.m - Compute the Hessian matrix for a packing.
%
% Author: Yuchen Xie
% Copyright (c) 2024-2026 Yuchen Xie.
%
% Harmonic / Hertz potential energy
%
% Input / output:
% simInfo:      current state; this function appends the Hessian matrix

%% Settings.
n = size(simInfo.xyz, 1);
if log2(n) ~= floor(log2(n))
    simInfo.xyz = [simInfo.xyz; [0, 0, 0]];
end

rcut = simInfo.radius * 2.0;
[N, dim] = size(simInfo.xyz);
simInfo.hessian = zeros([dim*N, dim*N]);
I3 = eye(dim);

%% Calculate Hessian.
dr = zeros([N-1, dim+2]);
for iatom = 1:N
    dr(:, 1) = simInfo.xyz(iatom, 1) - simInfo.xyz([1:iatom-1, iatom+1:N], 1);
    dr(:, 2) = simInfo.xyz(iatom, 2) - simInfo.xyz([1:iatom-1, iatom+1:N], 2);
    dr(:, 3) = simInfo.xyz(iatom, 3) - simInfo.xyz([1:iatom-1, iatom+1:N], 3);
    % periodic boundary condition
    dr(:, 1) = dr(:, 1) - round(dr(:, 1)/simInfo.simBoxLen(1)) * simInfo.simBoxLen(1);
    dr(:, 2) = dr(:, 2) - round(dr(:, 2)/simInfo.simBoxLen(2)) * simInfo.simBoxLen(2);
    dr(:, 3) = dr(:, 3) - round(dr(:, 3)/simInfo.simBoxLen(3)) * simInfo.simBoxLen(3);
    dr(:, 4) = sqrt(dr(:, 1).^2 + dr(:, 2).^2 + dr(:, 3).^2);
    dr(:, 5) = [1:iatom - 1, iatom + 1:N];

    cinfo = dr(dr(:, 4) <= rcut, :);
    if isempty(cinfo)
        continue;
    end

    cinfo(:, 1) = cinfo(:, 1)./cinfo(:, 4);
    cinfo(:, 2) = cinfo(:, 2)./cinfo(:, 4);
    cinfo(:, 3) = cinfo(:, 3)./cinfo(:, 4);

    rdivsig = cinfo(:, 4) / rcut;                           % d_{ij}/d

    % Harmonic
    tij = -(1.0 - rdivsig) / rcut;                          % -(1-d_{ij}/d)/d
    sij = 1.0/rcut/rcut * ones([size(cinfo, 1), 1]);        % 1/d^2

    % Hertz
    % tij = -(1.0 - rdivsig).^(3/2) / rcut;                         % -(1-d_{ij}/d)^(3/2) / d
    % sij = 1.5 * (1.0 - rdivsig).^(1/2) / rcut / rcut;             % 3 * (1-d_{ij}/d)^(1/2) / (2 * d^2)

    for jth = 1:size(cinfo, 1)
        jatom = cinfo(jth, 5);       % the number of the particle that is contacted with iatom
        
        Nij = cinfo(jth, 1:3)' * cinfo(jth, 1:3);
        hij = -(sij(jth) - tij(jth)/cinfo(jth, 4)) * Nij - tij(jth)/cinfo(jth, 4) * I3;

        simInfo.hessian(dim*(iatom - 1)+1:dim*iatom, dim*(jatom - 1)+1:dim*jatom) = hij;
        simInfo.hessian(dim*(iatom - 1)+1:dim*iatom, dim*(iatom - 1)+1:dim*iatom) = simInfo.hessian(dim*(iatom - 1)+1:dim*iatom, dim*(iatom - 1)+1:dim*iatom) - hij;
    end

end

%% Fix one particle.
simInfo.hessian = simInfo.hessian(1:3*N-3, 1:3*N-3);
simInfo.xyz = simInfo.xyz(1:N-1, :);

end
