% jampel_calcsubD.m - Compute pairwise distances between minima and build
%                     the corresponding sub-d matrix.
%
% Author: Yuchen Xie
% Copyright (c) 2024-2026 Yuchen Xie.
%
% This script:
% 1. Loads index-0 saddles from ../results/saddle_order_0.mat.
% 2. Computes all pairwise configuration distances.
% 3. Builds a weighted graph and its minimum spanning tree.
% 4. Saves the full distance matrix d_mat and the derived sub_d_mat.

addpath('../main');

%% Load minima data and allocate storage.
fprintf('Loading data...\n');
tic;
load('../results/saddle_order_0.mat');
N = size(saddle_order_0, 2);
weights = zeros(1, N*(N-1)/2);
s = zeros(1, N*(N-1)/2);
t = zeros(1, N*(N-1)/2);     % an edge from s to t
time = toc;
fprintf('Finished, time is %i \n', time);

%% Compute all pairwise distances between minima.
fprintf('Calculating distance \n');
tic;
parfor n = 1:N*(N-1)/2
    [i, j] = idx_trans(N, n);
    s(n) = i;
    t(n) = j;
    weights(n) = jampel_calcDistance(saddle_order_0(i), saddle_order_0(j));
end
time = toc;
fprintf('Finished, time is %i \n', time);

%% Build the complete weighted graph and its minimum spanning tree.
fprintf('minimum spanning tree \n');
tic;
G = graph(s, t, weights);
T = minspantree(G);
time = toc;
fprintf('Finished, time is %i \n', time);

%% Compute the all-pairs distance matrix and the sub-d matrix.
% sub_d_mat(i, j) is the largest edge weight along the shortest-tree path
% between minima i and j.
fprintf('Calculating sub_d \n');
tic;
d_mat = distances(G);
save('../results/d_mat.mat', 'd_mat');

sub_d_mat = zeros(N);
for i = 1:N-1
    for j = i+1:N
        P = shortestpath(T, i, j);
        nn = length(P);
        dd = zeros(1, nn-1);
        for k = 1:nn-1
            dd(k) = d_mat(P(k), P(k+1));
        end
        sub_d_mat(i, j) = max(dd);
    end
end
sub_d_mat = sub_d_mat + sub_d_mat';
time = toc;
fprintf('Finished, time is %i \n', time);

fprintf('Saving data...\n');
save('../results/d_mat.mat', 'd_mat');
save('../results/sub_d_mat.mat', 'sub_d_mat');
fprintf('The end. \n');

%% Convert a linear edge index to an upper-triangular matrix index pair.
function [i, j] = idx_trans(N, n)
i = 1;
while i < N
    if n <= N-i
        j = n + i;
        break
    else
        n = n - (N - i);
        i = i + 1;
    end
end
end
