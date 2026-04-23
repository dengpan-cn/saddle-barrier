% main.m - Starting from an index-4 saddle, enumerate connected minima
%          through repeated downward HiSD searches.
%
% Author: Yuchen Xie
% This workflow script was written by Yuchen Xie and uses project
% functions that include code adapted from earlier implementations by
% Jianyuan Yin.
% Copyright (c) 2024-2026 Yuchen Xie.
%
% This script loads previously generated saddle data from ../results and
% traverses the hierarchy from index 4 down to index 0 while recording
% child connections between adjacent saddle orders.

warning off

close all
clc

load('../results/saddle_order_4.mat');
load('../results/saddle_order_3.mat');
load('../results/saddle_order_2.mat');
load('../results/saddle_order_1.mat');
load('../results/saddle_order_0.mat');

l4 = length(saddle_order_4);

if l4 ~= 1
    error('More than one index-4 saddle point!!!');
end

%% Index-4 saddle to index-3 saddle.

for k = 1:1
    simInfo4 = saddle_order_4(k);
    s4 = jampel_calcForce(simInfo4);
    Index3 = [];
    Edge = [];
    for ind = 1:4           % stands for 4 unstable directions
        ii = 1;             % positive disturbance direction
        for i = 1:8         % search for 8 times
            fprintf('\n  4-3 new turn  \n');
            disp([k, ind, ii, i]);
            [simInfo3, iter] = down_search_basic(simInfo4, 4, 3, ind, ii);
            if iter == 100000 || iter == -1
                continue;
            else
                simInfo3 = jampel_cvtRun2ArXiv(simInfo3);
                s3 = jampel_calcForce(simInfo3);
                simInfo3 = jampel_calcHessian(simInfo3);
                lambda = eig(simInfo3.hessian);

                % make sure that the new configuration is exact index-3, and the energy is lower than the previous index-4 saddle
                if lambda(1) < -1e-8 && lambda(2) < -1e-8 && lambda(3) < -1e-8 && lambda(4) > 1e-8 && s4.ePair > s3.ePair
                    fprintf('yes\n');
                    simInfo3 = jampel_cvtRun2ArXiv(simInfo3);
                    idx_branch = Check_struct(Index3, simInfo3);
                    idx_all = Check_struct(saddle_order_3, simInfo3);
                    if idx_branch == 0
                        Index3 = [Index3, simInfo3];
                        if idx_all == 0
                            saddle_order_3 = [saddle_order_3, simInfo3];
                            idx_all = length(saddle_order_3);
                        else
                            fprintf('position in saddle_order_3 = %i \n', idx_all);
                        end
                        Edge = [Edge, idx_all];
                    else
                        fprintf('position in Index3 = %i \n', idx_branch);
                    end
                else
                    fprintf('no\n');
                    disp(lambda(1:6));
                    continue;
                end
            end
        end
        ii = -1;                % negative disturbance direction
        for i = 1:8             % search for 8 times
            fprintf('\n  4-3 new turn  \n');
            disp([k, ind, ii, i]);
            [simInfo3, iter] = down_search_basic(simInfo4, 4, 3, ind, ii);
            if iter == 100000 || iter == -1
                continue;
            else
                simInfo3 = jampel_cvtRun2ArXiv(simInfo3);
                s3 = jampel_calcForce(simInfo3);
                simInfo3 = jampel_calcHessian(simInfo3);
                lambda = eig(simInfo3.hessian);

                % make sure that the new configuration is exact index-3, and the energy is lower than the previous index-4 saddle
                if lambda(1) < -1e-8 && lambda(2) < -1e-8 && lambda(3) < -1e-8 && lambda(4) > 1e-8 && s4.ePair > s3.ePair
                    fprintf('yes\n');
                    simInfo3 = jampel_cvtRun2ArXiv(simInfo3);
                    idx_branch = Check_struct(Index3, simInfo3);
                    idx_all = Check_struct(saddle_order_3, simInfo3);
                    if idx_branch == 0
                        Index3 = [Index3, simInfo3];
                        if idx_all == 0
                            saddle_order_3 = [saddle_order_3, simInfo3];
                            idx_all = length(saddle_order_3);
                        else
                            fprintf('position in saddle_order_3 = %i \n', idx_all);
                        end
                        Edge = [Edge, idx_all];
                    else
                        fprintf('position in Index3 = %i \n', idx_branch);
                    end
                else
                    fprintf('no\n');
                    disp(lambda(1:6));
                    continue;
                end
            end
        end
    end
    saddle_order_4(k).childrenList = sort(Edge);
    
    save('../results/saddle_order_3.mat', 'saddle_order_3');
    save('../results/saddle_order_4.mat', 'saddle_order_4');
end

%% Save intermediate data.
save('../results/saddle_order_3.mat', 'saddle_order_3');
save('../results/saddle_order_4.mat', 'saddle_order_4');

%% Index-3 saddle to index-2 saddle.
l3 = length(saddle_order_3);

for k = 1:l3
    simInfo3 = saddle_order_3(k);
    s3 = jampel_calcForce(simInfo3);
    Index2 = [];
    Edge = [];
    for ind = 1:3           % stands for 3 unstable directions
        ii = 1;             % positive disturbance direction
        for i = 1:8         % search for 8 times
            fprintf('\n  3-2 new turn  \n');
            fprintf('%d / %d, %d, %d, %d \n', k, l3, ind, ii, i);
            [simInfo2, iter] = down_search_basic(simInfo3, 3, 2, ind, ii);
            if iter == 100000 || iter == -1
                continue;
            else
                simInfo2 = jampel_cvtRun2ArXiv(simInfo2);
                s2 = jampel_calcForce(simInfo2);
                simInfo2 = jampel_calcHessian(simInfo2);
                lambda = eig(simInfo2.hessian);

                % make sure that the new configuration is exact index-2, and the energy is lower than the previous index-3 saddle
                if lambda(1) < -1e-8 && lambda(2) < -1e-8 && lambda(3) > 1e-8 && s3.ePair > s2.ePair
                    fprintf('yes\n');
                    simInfo2 = jampel_cvtRun2ArXiv(simInfo2);
                    idx_branch = Check_struct(Index2, simInfo2);
                    idx_all = Check_struct(saddle_order_2, simInfo2);
                    if idx_branch == 0
                        Index2 = [Index2, simInfo2];
                        if idx_all == 0
                            saddle_order_2 = [saddle_order_2, simInfo2];
                            idx_all = length(saddle_order_2);
                        else
                            fprintf('position in saddle_order_2 = %i \n', idx_all);
                        end
                        Edge = [Edge, idx_all];
                    else
                        fprintf('position in Index2 = %i \n', idx_branch);
                    end
                else
                    fprintf('no\n');
                    disp(lambda(1:6));
                    continue;
                end
            end
        end
        ii = -1;            % negative disturbance direction
        for i = 1:8         % search for 8 times
            fprintf('\n  3-2 new turn  \n');
            fprintf('%d / %d, %d, %d, %d \n', k, l3, ind, ii, i);
            [simInfo2, iter] = down_search_basic(simInfo3, 3, 2, ind, ii);
            if iter == 100000 || iter == -1
                continue;
            else
                simInfo2 = jampel_cvtRun2ArXiv(simInfo2);
                s2 = jampel_calcForce(simInfo2);
                simInfo2 = jampel_calcHessian(simInfo2);
                lambda = eig(simInfo2.hessian);

                % make sure that the new configuration is exact index-2, and the energy is lower than the previous index-3 saddle
                if lambda(1) < -1e-8 && lambda(2) < -1e-8 && lambda(3) > 1e-8 && s3.ePair > s2.ePair
                    fprintf('yes\n');
                    simInfo2 = jampel_cvtRun2ArXiv(simInfo2);
                    idx_branch = Check_struct(Index2, simInfo2);
                    idx_all = Check_struct(saddle_order_2, simInfo2);
                    if idx_branch == 0
                        Index2 = [Index2, simInfo2];
                        if idx_all == 0
                            saddle_order_2 = [saddle_order_2, simInfo2];
                            idx_all = length(saddle_order_2);
                        else
                            fprintf('position in saddle_order_2 = %i \n', idx_all);
                        end
                        Edge = [Edge, idx_all];
                    else
                        fprintf('position in Index2 = %i \n', idx_branch);
                    end
                else
                    fprintf('no\n');
                    disp(lambda(1:6));
                    continue;
                end
            end
        end
    end
    saddle_order_3(k).childrenList = sort(Edge);
    
    save('../results/saddle_order_2.mat', 'saddle_order_2');
    save('../results/saddle_order_3.mat', 'saddle_order_3');
end

%% Save intermediate data.
save('../results/saddle_order_2.mat', 'saddle_order_2');
save('../results/saddle_order_3.mat', 'saddle_order_3');

%% Index-2 saddle to index-1 saddle.
l2 = length(saddle_order_2);

for k = 1:l2
    simInfo2 = saddle_order_2(k);
    s2 = jampel_calcForce(simInfo2);
    Index1 = [];
    Edge = [];
    for ind = 1:2           % stands for 2 unstable directions
        ii = 1;             % positive disturbance direction
        for i = 1:8         % search for 8 times
            fprintf('\n  2-1 new turn  \n');
            fprintf('%d / %d, %d, %d, %d \n', k, l2, ind, ii, i);
            [simInfo1, iter] = down_search_basic(simInfo2, 2, 1, ind, ii);
            if iter == 100000 || iter == -1
                continue;
            else
                simInfo1 = jampel_cvtRun2ArXiv(simInfo1);
                s1 = jampel_calcForce(simInfo1);
                simInfo1 = jampel_calcHessian(simInfo1);
                lambda = eig(simInfo1.hessian);

                % make sure that the new configuration is exact index-1, and the energy is lower than the previous index-2 saddle
                if lambda(1) < -1e-8 && lambda(2) > 1e-8 && s2.ePair > s1.ePair
                    fprintf('yes\n');
                    simInfo1 = jampel_cvtRun2ArXiv(simInfo1);
                    idx_branch = Check_struct(Index1, simInfo1);
                    idx_all = Check_struct(saddle_order_1, simInfo1);
                    if idx_branch == 0
                        Index1 = [Index1, simInfo1];
                        if idx_all == 0
                            saddle_order_1 = [saddle_order_1, simInfo1];
                            idx_all = length(saddle_order_1);
                        else
                            fprintf('position in saddle_order_1 = %i \n', idx_all);
                        end
                        Edge = [Edge, idx_all];
                    else
                        fprintf('position in Index1 = %i \n', idx_branch);
                    end
                else
                    fprintf('no\n');
                    disp(lambda(1:6));
                    continue;
                end
            end
        end
        ii = -1;            % negative disturbance direction
        for i = 1:8         % search for 8 times
            fprintf('\n  2-1 new turn  \n');
            fprintf('%d / %d, %d, %d, %d \n', k, l2, ind, ii, i);
            [simInfo1, iter] = down_search_basic(simInfo2, 2, 1, ind, ii);
            if iter == 100000 || iter == -1
                continue;
            else
                simInfo1 = jampel_cvtRun2ArXiv(simInfo1);
                s1 = jampel_calcForce(simInfo1);
                simInfo1 = jampel_calcHessian(simInfo1);
                lambda = eig(simInfo1.hessian);

                % make sure that the new configuration is exact index-1, and the energy is lower than the previous index-2 saddle
                if lambda(1) < -1e-8 && lambda(2) > 1e-8 && s2.ePair > s1.ePair
                    fprintf('yes\n');
                    simInfo1 = jampel_cvtRun2ArXiv(simInfo1);
                    idx_branch = Check_struct(Index1, simInfo1);
                    idx_all = Check_struct(saddle_order_1, simInfo1);
                    if idx_branch == 0
                        Index1 = [Index1, simInfo1];
                        if idx_all == 0
                            saddle_order_1 = [saddle_order_1, simInfo1];
                            idx_all = length(saddle_order_1);
                        else
                            fprintf('position in saddle_order_1 = %i \n', idx_all);
                        end
                        Edge = [Edge, idx_all];
                    else
                        fprintf('position in Index1 = %i \n', idx_branch);
                    end
                else
                    fprintf('no\n');
                    disp(lambda(1:6));
                    continue;
                end
            end
        end
    end
    saddle_order_2(k).childrenList = sort(Edge);
    
    save('../results/saddle_order_1.mat', 'saddle_order_1');
    save('../results/saddle_order_2.mat', 'saddle_order_2');
end

%% Save intermediate data.
save('../results/saddle_order_1.mat', 'saddle_order_1');
save('../results/saddle_order_2.mat', 'saddle_order_2');

%% Index-1 saddle to index-0 saddle.
l1 = length(saddle_order_1);

for k = 1:l1
    simInfo1 = saddle_order_1(k);
    s1 = jampel_calcForce(simInfo1);
    Index0 = [];
    Edge = [];
    % There is only one unstable direction.
    ii = 1;             % positive disturbance direction
    for i = 1:8         % search for 8 times
        fprintf('\n  1-0 new turn  \n');
        fprintf('%d / %d, %d, %d \n', k, l1, ii, i);
        [simInfo0, iter] = down_search_basic(simInfo1, 1, 0, 1, ii);
        if iter == 100000 || iter == -1
            continue;
        else
            simInfo0 = jampel_cvtRun2ArXiv(simInfo0);
            s0 = jampel_calcForce(simInfo0);
            simInfo0 = jampel_calcHessian(simInfo0);
            lambda = eig(simInfo0.hessian);

            % make sure that the new configuration is exact index-0, and the energy is lower than the previous index-1 saddle
            if lambda(1) > 1e-8 && s1.ePair > s0.ePair
                fprintf('yes\n');
                simInfo0 = jampel_cvtRun2ArXiv(simInfo0);
                idx_branch = Check_struct(Index0, simInfo0);
                idx_all = Check_struct(saddle_order_0, simInfo0);
                if idx_branch == 0
                    Index0 = [Index0, simInfo0];
                    if idx_all == 0
                        saddle_order_0 = [saddle_order_0, simInfo0];
                        idx_all = length(saddle_order_0);
                    else
                        fprintf('position in saddle_order_0 = %i \n', idx_all);
                    end
                    Edge = [Edge, idx_all];
                else
                    fprintf('position in Index0 = %i \n', idx_branch);
                end
            else
                fprintf('no\n');
                disp(lambda(1:6));
                continue;
            end
        end
    end
    ii = -1;            % negative disturbance direction
    for i = 1:8         % search for 8 times
        fprintf('\n  1-0 new turn  \n');
        fprintf('%d / %d, %d, %d \n', k, l1, ii, i);
        [simInfo0, iter] = down_search_basic(simInfo1, 1, 0, 1, ii);
        if iter == 100000 || iter == -1
            continue;
        else
            simInfo0 = jampel_cvtRun2ArXiv(simInfo0);
            s0 = jampel_calcForce(simInfo0);
            simInfo0 = jampel_calcHessian(simInfo0);
            lambda = eig(simInfo0.hessian);

            % make sure that the new configuration is exact index-0, and the energy is lower than the previous index-1 saddle
            if lambda(1) > 1e-8 && s1.ePair > s0.ePair
                fprintf('yes\n');
                simInfo0 = jampel_cvtRun2ArXiv(simInfo0);
                idx_branch = Check_struct(Index0, simInfo0);
                idx_all = Check_struct(saddle_order_0, simInfo0);
                if idx_branch == 0
                    Index0 = [Index0, simInfo0];
                    if idx_all == 0
                        saddle_order_0 = [saddle_order_0, simInfo0];
                        idx_all = length(saddle_order_0);
                    else
                        fprintf('position in saddle_order_0 = %i \n', idx_all);
                    end
                    Edge = [Edge, idx_all];
                else
                    fprintf('position in Index0 = %i \n', idx_branch);
                end
            else
                fprintf('no\n');
                disp(lambda(1:6));
                continue;
            end
        end
    end
    saddle_order_1(k).childrenList = sort(Edge);

    save('../results/saddle_order_0.mat', 'saddle_order_0');
    save('../results/saddle_order_1.mat', 'saddle_order_1');
end

%% Save final data.
save('../results/saddle_order_0.mat', 'saddle_order_0');
save('../results/saddle_order_1.mat', 'saddle_order_1');

fprintf('\n The end. \n');
