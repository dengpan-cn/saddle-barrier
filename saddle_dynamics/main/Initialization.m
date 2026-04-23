% Initialization.m - Generate a seed index-4 saddle by upsearch HiSD.
%
% Author: Yuchen Xie
% This workflow script was written by Yuchen Xie and uses project
% functions that include code adapted from earlier implementations by
% Jianyuan Yin.
% Copyright (c) 2024-2026 Yuchen Xie.
%
% This script:
% 1. Samples random initial packings.
% 2. Searches for a high-energy minimum.
% 3. Climbs successively to index-1 through index-4 saddles.
% 4. Saves the resulting saddle data under ../results.

%% Configuration settings.
phi = 0.655;        % fix a fraction
N = 64;             % fix an N (number of particles)
r = (3*phi/(4*pi*N))^(1/3);

%% Generate index-0 saddles (i.e. minima).
saddle_order_0(1, 100) = struct('xyz', [], 'simBoxLen', [], 'radius', [], ...
                                'childrenList', [], 'ratIdx', [], 'buckIdx', []);

idx = 1;
while idx <= 100
    disp(idx);
    
    simInfo = struct();
    simInfo.xyz = rand(N-1, 3) - 0.5;
    simInfo.simBoxLen = [1, 1, 1];
    simInfo.radius = r;

    [simInfo0, iter] = down_search_basic(simInfo, 1, 0, 1, 1);
    if ~(iter == -1 || iter == 100000)
        simInfo0 = jampel_cvtRun2ArXiv(simInfo0);
        saddle_order_0(idx) = simInfo0;
        idx = idx + 1;
    end
end

e = zeros(1, idx-1);
for i = 1:idx-1
    s = jampel_calcForce(saddle_order_0(i));
    e(i) = s.ePair;
end

[~, maxidx] = max(e);
simInfo0 = saddle_order_0(maxidx);      % Choose the minimum with the highest energy for upward search.
fprintf("found a minima\n\n");

%% Generate index-1 saddles.
saddle_order_1(1, 10)=struct('xyz', [], 'simBoxLen', [], 'radius', [], ...
                             'childrenList', [], 'ratIdx', [], 'buckIdx', []);

idx = 1;
while idx <= 10
    disp(idx);
    
    [simInfo1, iter] = up_search(simInfo0, 0, 1);
    if ~(iter == -1 || iter == 100000)
        simInfo1 = jampel_cvtRun2ArXiv(simInfo1);
        saddle_order_1(idx) = simInfo1;
        idx = idx + 1;
    end
end

e = zeros(1, idx-1);
for i = 1:idx-1
    s = jampel_calcForce(saddle_order_1(i));
    e(i) = s.ePair;
end

[~, maxidx] = max(e);
simInfo1 = saddle_order_1(maxidx);
fprintf("found a 1-index saddle point\n\n");

%% Generate index-2 saddles.
saddle_order_2(1, 10)=struct('xyz', [], 'simBoxLen', [], 'radius', [], ...
                             'childrenList', [], 'ratIdx', [], 'buckIdx', []);

idx = 1;
while idx <= 10
    disp(idx);
    
    [simInfo2, iter] = up_search(simInfo1, 1, 2);
    if ~(iter == -1 || iter == 100000)
        simInfo2 = jampel_cvtRun2ArXiv(simInfo2);
        saddle_order_2(idx) = simInfo2;
        idx = idx + 1;
    end
end

e = zeros(1, idx-1);
for i = 1:idx-1
    s = jampel_calcForce(saddle_order_2(i));
    e(i) = s.ePair;
end

[~, maxidx] = max(e);
simInfo2 = saddle_order_2(maxidx);
fprintf("found a 2-index saddle point\n\n");

%% Generate index-3 saddles.
saddle_order_3(1, 10)=struct('xyz', [], 'simBoxLen', [], 'radius', [], ...
                             'childrenList', [], 'ratIdx', [], 'buckIdx', []);

idx = 1;
while idx <= 10
    disp(idx);
    
    [simInfo3, iter] = up_search(simInfo2, 2, 3);
    if ~(iter == -1 || iter == 100000)
        simInfo3 = jampel_cvtRun2ArXiv(simInfo3);
        saddle_order_3(idx) = simInfo3;
        idx = idx + 1;
    end
end

e = zeros(1, idx-1);
for i = 1:idx-1
    s = jampel_calcForce(saddle_order_3(i));
    e(i) = s.ePair;
end

[~, maxidx] = max(e);
simInfo3 = saddle_order_3(maxidx);
fprintf("found a 3-index saddle point\n\n");

%% Generate index-4 saddles.
saddle_order_4(1, 10)=struct('xyz', [], 'simBoxLen', [], 'radius', [], ...
                             'childrenList', [], 'ratIdx', [], 'buckIdx', []);

idx = 1;
while idx <= 10
    disp(idx);
    
    [simInfo4, iter] = up_search(simInfo3, 3, 4);
    if ~(iter == -1 || iter == 100000)
        simInfo4 = jampel_cvtRun2ArXiv(simInfo4);
        saddle_order_4(idx) = simInfo4;
        idx = idx + 1;
    end
end

e = zeros(1, idx-1);
for i = 1:idx-1
    s = jampel_calcForce(saddle_order_4(i));
    e(i) = s.ePair;
end

[~, maxidx] = max(e);
simInfo4 = saddle_order_4(maxidx);
fprintf("found a 4-index saddle point\n\n");

%% Save data.
saddle_order_4 = simInfo4;
saddle_order_3 = [];
saddle_order_2 = [];
saddle_order_1 = [];
saddle_order_0 = [];

if ~exist('../results', 'dir')
    mkdir('../results');
end

save('../results/saddle_order_0.mat', 'saddle_order_0');
save('../results/saddle_order_1.mat', 'saddle_order_1');
save('../results/saddle_order_2.mat', 'saddle_order_2');
save('../results/saddle_order_3.mat', 'saddle_order_3');
save('../results/saddle_order_4.mat', 'saddle_order_4');

fprintf('\n The end. \n');
