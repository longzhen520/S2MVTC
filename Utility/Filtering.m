%Filtering
load('results.mat');

[bestACC, bestIdx] = max(results(:, 6));
bestParameters_ACC = results(bestIdx, :);

[bestNMI, bestIdx] = max(results(:, 7));
bestParameters_NMI = results(bestIdx, :);

[bestPurity, bestIdx] = max(results(:, 8));
bestParameters_Purity = results(bestIdx, :);

[bestTime, bestIdx] = min(results(:, 9));
bestParameters_Time = results(bestIdx, :);


figure;
plot(results(:, 6));
title('ACC under different parameter combinations');
xlabel('Parameter Combination Index');
ylabel('ACC');

figure;
plot(results(:, 7));
title('NMI under different parameter combinations');
xlabel('Parameter Combination Index');
ylabel('NMI');

figure;
plot(results(:, 8));
title('Purity under different parameter combinations');
xlabel('Parameter Combination Index');
ylabel('Purity');

figure;
plot(results(:, 9));
title('Time under different parameter combinations');
xlabel('Parameter Combination Index');
ylabel('Time');



% 例如，假设你为每个指标分配了相等的权重
weights = [0.25, 0.25, 0.25, 0.25];  % ACC, NMI, Purity, Time

% 归一化时间（假设时间越短越好）
normalized_time = (max(results(:, 9)) - results(:, 9)) / (max(results(:, 9)) - min(results(:, 9)));

% 计算总分
scores = weights(1)*results(:, 6) + weights(2)*results(:, 7) + weights(3)*results(:, 8) + weights(4)*normalized_time;

% 找到最高分
[bestScore, bestIdx] = max(scores);
bestParameters = results(bestIdx, :);


paretoSet = findParetoFront([results(:, 6:8), -results(:, 9)]);  % 优化Time为负数，因为时间越短越好参数
paretoResults = results(paretoSet, :);


% 输出每个指标的最优参数组合
fprintf('Best Parameters for ACC: Beta = %.5f, Gamma = %.5f, Lambda = %.5f, Alpha_1 = %.5f, Alpha_2 = %.5f, ACC = %.5f, NMI = %.5f, Purity = %.5f, Time = %.5f\n', bestParameters_ACC);
fprintf('Best Parameters for NMI: Beta = %.5f, Gamma = %.5f, Lambda = %.5f, Alpha_1 = %.5f, Alpha_2 = %.5f, ACC = %.5f, NMI = %.5f, Purity = %.5f, Time = %.5f\n', bestParameters_NMI);
fprintf('Best Parameters for Purity: Beta = %.5f, Gamma = %.5f, Lambda = %.5f, Alpha_1 = %.5f, Alpha_2 = %.5f, ACC = %.5f, NMI = %.5f, Purity = %.5f, Time = %.5f\n', bestParameters_Purity);
fprintf('Best Parameters for Time: Beta = %.5f, Gamma = %.5f, Lambda = %.5f, Alpha_1 = %.5f, Alpha_2 = %.5f, ACC = %.5f, NMI = %.5f, Purity = %.5f, Time = %.5f\n', bestParameters_Time);

% 输出加权评分的最优参数组合
fprintf('Best Parameters with Weighted Score: Beta = %.5f, Gamma = %.5f, Lambda = %.5f, Alpha_1 = %.5f, Alpha_2 = %.5f, ACC = %.5f, NMI = %.5f, Purity = %.5f, Time = %.5f\n', bestParameters);

% 输出帕累托最优解
for i = 1:size(paretoResults, 1)
    fprintf('Pareto Solution %d: Beta = %.5f, Gamma = %.5f, Lambda = %.5f, Alpha_1 = %.5f, Alpha_2 = %.5f, ACC = %.5f, NMI = %.5f, Purity = %.5f, Time = %.5f\n', i, paretoResults(i, :));
end



% 加载数据
load('results.mat');

% 获取第 924 行的参数
params_924 = results(924, :);

% 输出第 924 行的参数
fprintf('Parameters on row 924: Beta = %.6f, Gamma = %.6f, Lambda = %.6f, Alpha_1 = %.6f, Alpha_2 = %.6f, ACC = %.6f, NMI = %.6f, Purity = %.6f, Time = %.6f\n', params_924);


