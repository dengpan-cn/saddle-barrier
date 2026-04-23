% jampel_calcsubE.m - Compute barrier-based connectivity between minima and
%                     build the corresponding sub-E matrix.
%
% Author: Yuchen Xie
% Copyright (c) 2024-2026 Yuchen Xie.
%
% This script:
% 1. Loads archived index-0 and index-1 saddles from ../results.
% 2. Builds a graph whose edges are induced by shared index-1 saddles.
% 3. Restricts the analysis to the largest connected component ("hub").
% 4. Computes sub_E_mat for the hub and saves the restricted sub_d matrix.

addpath('../main');

%% Load archived minima and index-1 saddles.
load('../results/saddle_order_0.mat');
load('../results/saddle_order_1.mat');

n0 = length(saddle_order_0);
n1 = length(saddle_order_1);
count = 0;
N = n0*(n0-1) / 2;
s = zeros(1, N); t = zeros(1, N); weight = zeros(1, N);

%% Precompute the energy of each minimum.
E0 = zeros(1, n0);
for i = 1:n0
    s0 = jampel_calcForce(saddle_order_0(i));
    E0(i) = s0.ePair;
end

%% Calculate the energy barrier on a path from index-1 to index-0.
for i = 1:n1
    chList = saddle_order_1(i).childrenList;
    nch = length(chList);
    S1_temp = jampel_calcForce(saddle_order_1(i));
    Esaddle = S1_temp.ePair;

    for k = 1:nch
        for j = k+1:nch
            count = count + 1;
            s(count) = chList(k);
            t(count) = chList(j);
            weight(count) = Esaddle;
        end
    end
end
s = s(1:count);
t = t(1:count);
weight = weight(1:count);
H = graph(s, t, weight);

G = simplify(H, 'min');

%% Find the largest connected component in the barrier graph.
bins = conncomp(G);
bins_num = max(bins);
bn = zeros(1, bins_num);
for i = 1:bins_num
    bn(i) = length(find(bins == i));
end
[~, index] = max(bn);

hub_idx = find(bins == index);
hub_num = length(hub_idx);
fprintf('The number of the maximum hub is %i \n', hub_num);

%% Compute the sub-E matrix on the largest connected component.
% sub_E_mat(i, j) stores the largest barrier along the tree path between
% two minima, measured relative to the endpoint energies.
G = subgraph(G, hub_idx);
E_mat = distances(G);
T = minspantree(G);

sub_E_mat = zeros(hub_num);

for i = 1:hub_num-1
    for j = i+1:hub_num
        P = shortestpath(T, i, j);
        nn = length(P);
        dd = zeros(1, nn-1);
        for k = 1:nn-1
            dd(k) = E_mat(P(k), P(k+1));
        end
        md = max(dd);
        sub_E_mat(i, j) = max([md-E0(hub_idx(i)), md-E0(hub_idx(j))]);
    end
end
sub_E_mat = sub_E_mat + sub_E_mat';

%% Restrict the previously computed sub_d matrix to the same hub.
load('../results/sub_d_mat.mat');
sub_d_hub = sub_d_mat(hub_idx, hub_idx);
save('../results/sub_d_hub.mat', 'sub_d_hub');
save('../results/sub_E_mat.mat', 'sub_E_mat');

fprintf('The end. \n');
