function paretoSet = findParetoFront(data)
    % data: 每行是一个解，每列是一个目标。
    % paretoSet: 一个逻辑向量，标记哪些行是帕累托最优解。

    numSolutions = size(data, 1);
    paretoSet = true(numSolutions, 1);  % 初始假设所有解都是帕累托最优解

    for i = 1:numSolutions
        if paretoSet(i)  % 只检查仍然被认为是帕累托最优解的解
            for j = 1:numSolutions
                if i ~= j && all(data(j, :) <= data(i, :)) && any(data(j, :) < data(i, :))
                    paretoSet(i) = false;  % 找到了一个更好的解，所以解 i 不是帕累托最优解
                    break;
                end
            end
        end
    end
end

