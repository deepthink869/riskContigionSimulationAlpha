%% Financial Group Risk Contagion Simulation (t-Copula + Network Model)
% This script simulates default contagion within a corporate group via Monte Carlo,
% incorporating both copula-based default correlation and network-based contagion.
% It supports three input modes: built-in example data, interactive manual input, 
% or reading from an Excel file. The output includes:
% - Average number of initial defaults
% - Average number of total infected (defaulted) companies 
% - Probability of a large cascade (contagion affecting >= threshold companies)
% - Frequency each company triggers a large cascade
% - Estimated contagion probability matrix (empirical frequencies from simulation)
% Results can be displayed and optionally saved to an Excel file.

%% 1. Data Input
disp('请选择数据输入方式:');
disp('1 = 使用内置示例数据');
disp('2 = 手动输入数据 (提供公司属性和关系, 自动计算参数)');
disp('3 = 从Excel文件读取数据');
choice = input('请输入选项 [1/2/3]: ');
if isempty(choice)
    choice = 1;  % default to example data if no input
end

% Initialize variables
n = [];               % number of companies
defaultProb = [];     % default probability vector (annual PD for each company)
W = [];               % contagion adjacency matrix (n x n)
R = [];               % correlation matrix for copula (n x n)
df = [];              % t-Copula degrees of freedom

if choice == 1
    % === Built-in Example Data (10 companies) ===
    n = 10;
    % Default probabilities vector (based on assumed credit ratings):
    defaultProb = [0.005; 0.010; 0.020; 0.001; 0.015; 0.020; 0.003; 0.012; 0.008; 0.005];
    % Contagion network matrix W (10x10). We specify a few contagion links with probabilities:
    W = zeros(n);
    % Example assumptions for W:
    W(2,1) = 0.3;   % Company2 default -> Company1 defaults with 0.3 probability (e.g. 1 guaranteed 2's loans)
    W(3,2) = 0.2;   % Company3 default -> Company2 with 0.2 probability
    W(4,2) = 0.1;   % Company4 default -> Company2 with 0.1 probability (some minor link)
    W(3,5) = 0.5;   % Company3 default -> Company5 with 0.5 probability (e.g. 5 is a key supplier to 3 and loses payments)
    W(7,5) = 0.25;  % Company7 default -> Company5 with 0.25 probability
    W(8,6) = 0.15;  % Company8 default -> Company6 with 0.15 probability
    W(9,7) = 0.30;  % Company9 default -> Company7 with 0.30 probability
    W(10,7)= 0.40;  % Company10 default -> Company7 with 0.40 probability
    % Copula correlation matrix R (10x10):
    R = eye(n);
    R(1,2) = 0.5; R(2,1) = 0.5;   % Companies 1 and 2 have high correlation (e.g. same industry or group)
    R(3,4) = 0.4; R(4,3) = 0.4;   % Companies 3 and 4 moderately correlated 
    % (Other off-diagonals remain 0 or low implying minimal correlation)
    % t-Copula degrees of freedom:
    df = 5;
    fprintf('已载入内置示例参数：%d 家企业。\n', n);

elseif choice == 2
    % === Interactive Manual Input Mode ===
    n = input('请输入企业数量 n: ');
    if isempty(n) || n<=0
        error('企业数量必须是正整数.');
    end

    % Prompt for advanced input or direct input:
    disp('请选择违约概率(PD)输入方式:');
    disp('1 = 直接输入每家企业的违约概率');
    disp('2 = 输入信用评级或风险等级, 由系统映射为违约概率');
    pdMode = input('选项 [1/2]: ');
    if isempty(pdMode)
        pdMode = 2;  % default to rating-based input
    end

    defaultProb = zeros(n,1);
    if pdMode == 1
        % Direct PD input for each company
        fprintf('请输入每家企业的年违约概率 (0~1之间):\n');
        for i = 1:n
            p = input(sprintf(' 公司 %d 的违约概率: ', i));
            if isempty(p) || p < 0 || p > 1
                error('输入错误：公司 %d 的违约概率无效。请输0到1之间的数。', i);
            end
            defaultProb(i) = p;
        end
    else
        % Rating or qualitative input for each company -> map to PD
        % Define a mapping from ratings to PD (these are rough typical values):
        ratingPDMap = struct('AAA',0.0001,'AA',0.0005,'A',0.001,'BBB',0.005,'BB',0.02,'B',0.05,'CCC',0.15);
        % (If a rating is not in this map, we will prompt for a numeric PD.)
        fprintf('请输入每家企业的信用评级或风险等级 (AAA/AA/A/BBB/BB/B/CCC 等，或直接输入概率):\n');
        for i = 1:n
            resp = input(sprintf(' 公司 %d: ', i),'s');
            if isempty(resp)
                error('输入错误：未提供公司 %d 的评级或概率。', i);
            end
            resp = strtrim(resp);
            % Check if numeric input (the user can input a number as string)
            numVal = str2double(resp);
            if ~isnan(numVal)
                % If a numeric value was entered as string, use it directly
                if numVal < 0 || numVal > 1
                    error('公司 %d 的违约概率数值无效，请输入0~1之间。', i);
                end
                defaultProb(i) = numVal;
            else
                % If not purely numeric, interpret as rating code
                upperCode = upper(resp);
                if isfield(ratingPDMap, upperCode)
                    defaultProb(i) = ratingPDMap.(upperCode);
                else
                    % Unrecognized rating string: prompt for numeric probability
                    warning('未知的评级 "%s". 请直接提供公司 %d 的违约概率数值:', resp, i);
                    p = input(sprintf(' 公司 %d 的违约概率: ', i));
                    if isempty(p) || p < 0 || p > 1
                        error('输入错误：公司 %d 的违约概率无效。', i);
                    end
                    defaultProb(i) = p;
                end
            end
        end
        fprintf('已根据评级输入计算违约概率向量。\n');
    end

    % Input contagion relationships to build W matrix
    W = zeros(n);
    % 2.1 Guarantee (Credit) relationships
    numGuar = input('请输入担保关系的数量 (企业间互保/担保数): ');
    if isempty(numGuar)
        numGuar = 0;
    end
    if numGuar > 0
        disp('请逐条输入担保关系：[被担保企业 索赔担保企业 担保暴露比例]');
        disp('  例如输入 "2 1 0.3" 表示: 公司2违约 -> 公司1有0.3概率违约 (1为2的担保方,暴露30%)');
    end
    for k = 1:numGuar
        rel = input(sprintf(' 担保关系 %d: ', k),'s');
        if isempty(rel)
            error('输入错误：第 %d 条担保关系未提供。', k);
        end
        vals = sscanf(rel,'%f');
        if numel(vals) < 2
            error('输入格式错误：请提供 "i j [权重]" 三个值。');
        end
        i_idx = round(vals(1));  % guaranteed (borrower) company
        j_idx = round(vals(2));  % guarantor company
        if i_idx<1 || i_idx>n || j_idx<1 || j_idx>n
            error('担保关系索引超出范围: 公司 %d 或 %d 不存在。', i_idx, j_idx);
        end
        % Determine weight
        w_ij = 0.3;  % default weight if not provided
        if numel(vals) >= 3
            w_ij = vals(3);
            if w_ij <= 0
                warning('担保概率权重应为正值，已默认使用0.3。');
                w_ij = 0.3;
            end
            if w_ij > 1
                % If user gave a percentage >1 (like 50 meaning 50%), convert to fraction
                w_ij = w_ij/100;
            end
            if w_ij > 1
                w_ij = 1;  % cap at 100%
            end
        end
        W(i_idx, j_idx) = w_ij;
    end

    % 2.2 Supply chain relationships
    numSupply = input('请输入供应链传染关系的数量 (关键客户-供应商数): ');
    if isempty(numSupply)
        numSupply = 0;
    end
    if numSupply > 0
        disp('请逐条输入供应链关系：[客户 企业 供应商 企业 依赖比例]');
        disp('  例如输入 "3 5 0.5" 表示: 公司3违约 -> 公司5有0.5概率违约 (5有50%业务依赖于3)');
    end
    for k = 1:numSupply
        rel = input(sprintf(' 供应链关系 %d: ', k),'s');
        if isempty(rel)
            error('输入错误：第 %d 条供应链关系未提供。', k);
        end
        vals = sscanf(rel,'%f');
        if numel(vals) < 2
            error('输入格式错误：请提供 "i j [权重]" 三个值。');
        end
        i_idx = round(vals(1));  % customer (debtor) company
        j_idx = round(vals(2));  % supplier (creditor) company
        if i_idx<1 || i_idx>n || j_idx<1 || j_idx>n
            error('供应链关系索引超出范围: 公司 %d 或 %d 不存在。', i_idx, j_idx);
        end
        w_ij = 0.2;  % default weight if not provided
        if numel(vals) >= 3
            w_ij = vals(3);
            if w_ij <= 0
                warning('传染概率应为正值，已默认使用0.2。');
                w_ij = 0.2;
            end
            if w_ij > 1
                % If user gave a percentage >1, interpret as percentage
                w_ij = w_ij/100;
            end
            if w_ij > 1
                w_ij = 1;
            end
        end
        W(i_idx, j_idx) = w_ij;
    end

    % (Optional) 2.3 Other relationships (e.g. equity/group) could be added similarly...
    % For simplicity, we'll assume major relationships are covered by guarantee or supply chain.
    % Users can adjust W manually here if needed for any additional links.

    % 3. Copula correlation matrix R based on industry/group
    R = eye(n);
    % Get industry or group labels for each company:
    industryLabels = cell(n,1);
    disp('请输入每家企业的行业或所属组别标识:');
    for i = 1:n
        lbl = input(sprintf(' 公司 %d 行业/类别: ', i),'s');
        if isempty(lbl)
            lbl = sprintf('Company%d', i);  % if not provided, use unique default to avoid false grouping
        end
        industryLabels{i} = strtrim(lbl);
    end
    % Set baseline low correlation for all pairs:
    baseCorr = 0.1;
    R(1:n, 1:n) = baseCorr;  % set entire matrix to 0.1 (will fix diagonal and specific pairs next)
    for i = 1:n
        R(i,i) = 1;  % self-correlation
    end
    % Increase correlation for related companies:
    for i = 1:n
        for j = i+1:n
            % Same industry/group label -> high correlation
            if ~isempty(industryLabels{i}) && strcmp(industryLabels{i}, industryLabels{j})
                R(i,j) = 0.5;
                R(j,i) = 0.5;
            end
            % Direct financial link (from W) -> moderate correlation if not same industry
            if W(i,j) > 0 || W(j,i) > 0
                % If they already have higher correlation from industry, we keep that;
                % otherwise, set a moderate correlation for business connection:
                if R(i,j) < 0.3
                    R(i,j) = 0.3;
                    R(j,i) = 0.3;
                end
            end
        end
    end

    % t-Copula degrees of freedom
    df = input('请输入 t-Copula 的自由度参数 df (默认5): ');
    if isempty(df)
        df = 5;
    end
    fprintf('参数输入完毕，共 %d 家企业。\n', n);

elseif choice == 3
    % ==== 从 Excel 文件读取模式（跨平台稳定：readmatrix） ====
    fileName = input('请输入Excel文件名(含扩展名，例如 data.xlsx): ', 's');

    % 支持当前目录/桌面两地查找
    if exist(fileName, 'file') ~= 2
        altPath = fullfile(getenv('HOME'), 'Desktop', fileName);
        if exist(altPath, 'file') == 2
            fileName = altPath;
            disp(['已从桌面读取文件：', fileName]);
        else
            error('找不到文件: %s。请确认文件在当前目录或桌面。', fileName);
        end
    end

    % === 用 readmatrix 读取 3 个工作表 ===
    % Sheet1: PD 向量（N×1 或 1×N，纯数字，小数形式）
    defaultProb = readmatrix(fileName, 'Sheet', 1);   % 兼容 macOS/Windows/Linux
    defaultProb = defaultProb(:);                      % 转成列向量
    n = length(defaultProb);
    fprintf('从Excel读取到 %d 家企业的违约概率。\n', n);

    % Sheet2: W (N×N 传染概率矩阵)
    W = readmatrix(fileName, 'Sheet', 2);
    if ~ismatrix(W) || any(size(W) ~= [n, n])
        error('Excel中第二张表的 W 矩阵尺寸与企业数不符（应为 %d×%d）。', n, n);
    end

    % Sheet3: R (N×N Copula 相关矩阵)
    R = readmatrix(fileName, 'Sheet', 3);
    if ~ismatrix(R) || any(size(R) ~= [n, n])
        error('Excel中第三张表的 R 矩阵尺寸与企业数不符（应为 %d×%d）。', n, n);
    end

    % === 基本合法性检查（给出提示并自动修正一些小问题） ===
    if any(~isfinite(defaultProb))
        error('Sheet1（违约概率）存在非数值/NaN/Inf，请检查。');
    end
    if any(defaultProb < 0 | defaultProb > 1)
        error('Sheet1（违约概率）应在 [0,1] 之间（小数形式）。');
    end

    if any(~isfinite(W), 'all')
        error('Sheet2（W）存在非数值/NaN/Inf，请检查。');
    end
    if any(W(:) < 0 | W(:) > 1)
        error('Sheet2（W）元素必须在 [0,1] 之间。');
    end

    if any(~isfinite(R), 'all')
        error('Sheet3（R）存在非数值/NaN/Inf，请检查。');
    end
    % 若不完全对称或对角不为1，做温和修正并提示
    if max(abs(R - R.'), [], 'all') > 1e-8
        warning('R 非完全对称，已自动对称化处理。');
        R = (R + R.') / 2;
    end
    if max(abs(diag(R) - 1)) > 1e-8
        warning('R 对角线非 1，已自动修正为 1。');
        R(1:n+1:end) = 1;
    end

    % t-Copula 自由度
    df = input('请输入 t-Copula 的自由度参数 df: ');
    if isempty(df)
        df = 5;
    end
    fprintf('已从Excel读取参数，共 %d 家企业。\n', n);

% =======================================

    % Read defaultProb from sheet1:
    defaultProb = xlsread(fileName, 1);
    defaultProb = defaultProb(:);  % ensure column vector
    n = length(defaultProb);
    fprintf('从Excel读取到 %d 家企业的违约概率。\n', n);
    % Read W matrix from sheet2:
    W = xlsread(fileName, 2);
    if size(W,1) ~= n || size(W,2) ~= n
        error('Excel中第二张表的 W 矩阵尺寸与企业数不符 (%d×%d)。', size(W,1), size(W,2));
    end
    % Read R matrix from sheet3:
    R = xlsread(fileName, 3);
    if size(R,1) ~= n || size(R,2) ~= n
        error('Excel中第三张表的 R 矩阵尺寸与企业数不符 (%d×%d)。', size(R,1), size(R,2));
    end
    % Prompt for df parameter:
    df = input('请输入 t-Copula 的自由度参数 df: ');
    if isempty(df)
        df = 5;
    end
    fprintf('已从Excel读取参数，共 %d 家企业。\n', n);

else
    error('无效的选择，程序终止。');
end

%% 2. Simulation settings
numSim = 10000;    % number of Monte Carlo simulation runs (adjust as needed for precision)
rng(1);            % fix random seed for reproducibility

% Preallocate results storage
results_initial = zeros(numSim, 1);   % initial default count in each simulation
results_cascade = zeros(numSim, 1);   % final total infected count in each simulation
transmission_count = zeros(n, n);     % count of contagion occurrences along each edge
init_default_matrix = false(numSim, n); % record which firms default initially in each run (logical matrix)

%% 3. Monte Carlo Simulation
for sim = 1:numSim
    % Step 1: Generate one random shock vector from t-Copula (length n)
    U = copularnd('t', R, df, 1);  % 1xN random vector from t-copula with correlation R and df
    % Determine which firms default initially (if U_i < PD_i)
    init_default = (U < defaultProb');
    init_default_matrix(sim, :) = init_default;
    results_initial(sim) = sum(init_default);
    % Step 2: Contagion propagation over the network W
    infected = init_default;   % start with initially defaulted firms
    new_infected = true;
    while new_infected
        new_infected = false;
        % Spread from each currently infected node to its neighbors
        for i = 1:n
            if infected(i)
                for j = 1:n
                    if ~infected(j) && W(i,j) > 0
                        % Company i is defaulted and has a link to j with probability W(i,j)
                        if rand() < W(i,j)
                            infected(j) = true;
                            new_infected = true;
                            transmission_count(i,j) = transmission_count(i,j) + 1;
                        end
                    end
                end
            end
        end
    end
    % Step 3: Record final infected count for this run
    results_cascade(sim) = sum(infected);
end

% Compute empirical transmission probability matrix T from simulations
T = transmission_count / numSim;

%% 4. Result Summary
avg_initial = mean(results_initial);
avg_infected = mean(results_cascade);
max_infected = max(results_cascade);
min_infected = min(results_cascade);
fprintf('\n模拟完成。主要结果如下：\n');
disp(['平均初始违约企业数: ', num2str(avg_initial)]);
disp(['平均受感染企业数: ', num2str(avg_infected)]);
disp(['最大全部感染企业数: ', num2str(max_infected)]);
disp(['最小感染企业数: ', num2str(min_infected)]);

% Probability of a large contagion event (e.g., >= 5 companies infected)
threshold = 5;
large_outbreak_prob = mean(results_cascade >= threshold);
disp(['大规模连锁传染(>=', num2str(threshold),'家)发生概率: ', num2str(large_outbreak_prob)]);

% Frequency each company triggers a large cascade (as initial default)
large_trigger = zeros(1, n);
for i = 1:n
    count_i = 0;
    for sim = 1:numSim
        if init_default_matrix(sim, i) && results_cascade(sim) >= threshold
            count_i = count_i + 1;
        end
    end
    large_trigger(i) = count_i / numSim;
end
disp(['各节点触发大规模传染的频率 (按节点1-', num2str(n), '):']);
disp(num2str(large_trigger));

% Contagion probability matrix T
disp('传染概率矩阵 T (每条边发生传染的估计概率):');
disp(num2str(T));

%% 5. Save Results (optional)
saveFile = input('如需将结果保存到Excel, 请输入文件名(如 output.xlsx), 或直接回车跳过: ', 's');
if ~isempty(saveFile)
    % Ensure .xlsx extension
    if isempty(regexpi(saveFile, '\.xlsx$'))
        saveFile = [saveFile, '.xlsx'];
    end
    % If no path given, default to Desktop
    if ~contains(saveFile, filesep)
        homeDir = getenv('HOME');
        if isempty(homeDir)
            homeDir = pwd;
        end
        saveFile = fullfile(homeDir, 'Desktop', saveFile);
    end
    if exist(saveFile, 'file')
        delete(saveFile); % overwrite if exists
    end
    try
        % Summary sheet
        summaryData = {
            '平均初始违约数', avg_initial;
            '平均受感染企业数', avg_infected;
            '最大全部感染数', max_infected;
            '最小感染数', min_infected;
            ['>=', num2str(threshold), '家连锁违约概率'], large_outbreak_prob
        };
        writecell(summaryData, saveFile, 'Sheet', 'Summary', 'Range', 'A1');
        writecell({'各节点触发大规模传染频率:'}, saveFile, 'Sheet', 'Summary', 'Range', 'A7');
        writematrix(large_trigger, saveFile, 'Sheet', 'Summary', 'Range', 'A8');
        % TransmissionMatrix sheet
        T_cell = num2cell(T);
        T_with_labels = cell(n+1, n+1);
        T_with_labels{1,1} = '节点';
        for j = 1:n
            T_with_labels{1, j+1} = ['节点', num2str(j)];
            T_with_labels{j+1, 1} = ['节点', num2str(j)];
            for k = 1:n
                T_with_labels{j+1, k+1} = T_cell{j, k};
            end
        end
        writecell(T_with_labels, saveFile, 'Sheet', 'TransmissionMatrix', 'Range', 'A1');
        disp(['结果已保存至 Excel 文件：', saveFile]);
    catch ME
        disp('写入 Excel 时发生错误:');
        disp(ME.message);
    end
end

