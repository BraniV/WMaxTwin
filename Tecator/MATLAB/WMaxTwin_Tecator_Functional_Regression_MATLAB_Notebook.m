%% WMaxTwin for Tecator Functional Regression
% MATLAB notebook-style script for Example III in the WMaxTwin paper.
%
% This file is intentionally written as a sectioned MATLAB script using %%
% cells.  It can be opened directly in MATLAB, run cell by cell, and saved as
% a MATLAB Live Script (.mlx) if desired.
%
% The example compares split geometries for scalar-on-function regression on
% the Tecator near-infrared spectroscopy data.  The response is fat content.
% The split is over spectra, not over wavelengths.  Each spectrum keeps its
% internal grid, and Haar wavelet coefficients are computed within each
% spectrum.  WMaxTwin modifies the validation geometry by adding weighted
% wavelet-scale information to ordinary raw-spectrum MaxTwin.
%
% The workflow is:
%   1. Load StatLib Tecator spectra.
%   2. Interpolate 100 channels to 128 grid points.
%   3. Use samples 1:172 as calibration and 173:215 as external test.
%   4. Standardize spectra using calibration statistics only.
%   5. Compute an orthonormal Haar representation of every spectrum.
%   6. Construct raw distance D0 and wavelet distance Dw.
%   7. Construct nested WMaxTwin distances D_gamma = (1-gamma)D0 + gamma Dw.
%   8. Compare Random, SPlit, Twinning, DUPLEX, MaxTwin, and WMaxTwin splits.
%   9. Select jmax and ridge penalty by validation inside calibration.
%  10. Refit on all calibration samples and evaluate on the external test set.
%
% The code avoids toolboxes: Haar, ridge regression, split construction,
% distance calculations, CSV writing, and plotting are implemented directly.
% MATLAB's boxplot function is not used because it requires the Statistics and
% Machine Learning Toolbox in some installations.
%
% Expected output directory:
%     tecator_wmaxtwin_outputs_MATLAB
%
% Default run below uses MATLAB_REPS = 8 for speed.  Increase to 50 to match
% the paper-scale Monte Carlo more closely.

clear; close all; clc;

%% Configuration
DATA_URL = 'https://lib.stat.cmu.edu/datasets/tecator';
OUTDIR = fullfile(pwd, 'tecator_wmaxtwin_outputs_MATLAB');
DATA_PATH = fullfile(pwd, 'tecator_statlib.txt');
if ~exist(OUTDIR, 'dir'), mkdir(OUTDIR); end

MATLAB_REPS = 8;          % increase to 50 for paper-scale run
SEED = 20260619;
M_VAL = 43;               % validation size inside 172 calibration samples
ALPHAS = logspace(-4, 4, 13);
GAMMA_PATH = [0.00, 0.10, 0.25, 0.50, 0.75];

%% Load the Tecator data
% The package includes tecator_statlib.txt.  If it is missing, MATLAB will try
% to download the file from StatLib.

ensureTecatorData(DATA_PATH, DATA_URL);
[X, endpoints, wavelengths] = loadTecator(DATA_PATH);

% StatLib endpoints are moisture, fat, protein.  We use fat percentage.
y_fat = endpoints(:,2);

[X128, wl128] = interpToPower2(X, wavelengths, 128);

% Original protocol used in the paper: calibration/test split only.
% Samples 216:240 are omitted as extrapolation examples.
X_cal_raw  = X128(1:172, :);
y_cal      = y_fat(1:172);
X_test_raw = X128(173:215, :);
y_test     = y_fat(173:215);

fprintf('Tecator data loaded: %d calibration spectra, %d external test spectra.\n', ...
    size(X_cal_raw,1), size(X_test_raw,1));

%% Standardize spectra and compute Haar coefficients
% Standardization uses only the calibration set.  The same center and scale are
% applied to the external test spectra.

sc_curve = standardizeFit(X_cal_raw);
X_cal = standardizeApply(X_cal_raw, sc_curve);
X_test = standardizeApply(X_test_raw, sc_curve);

[H, blocks] = haarTransformMatrix(128);
C_cal = X_cal * H';
C_test = X_test * H';
JmaxGrid = 1:round(log2(size(X_cal,2)));

%% Construct raw MaxTwin and WMaxTwin distances
% D0 is the raw standardized-spectra geometry.  Dw is the wavelet geometry.
% The wavelet scale weights are computed on calibration data only.  They are
% not optimized using the external test set, and they are not raw energy weights.

D0 = normalizeDistance(pairwiseSqDists(X_cal));
weights = scaleRelevanceWeights(C_cal, y_cal, blocks);
Dw = waveletWeightedDistance(C_cal, blocks, weights);

equal_weights = weights;
for j = 1:numel(equal_weights)
    equal_weights(j).weight = 1 / numel(equal_weights);
end
Dw_equal = waveletWeightedDistance(C_cal, blocks, equal_weights);

D_gamma = cell(numel(GAMMA_PATH), 1);
for ii = 1:numel(GAMMA_PATH)
    g = GAMMA_PATH(ii);
    D_gamma{ii} = normalizeDistance((1 - g) * D0 + g * Dw);
end
D_gamma_equal = normalizeDistance(0.75 * D0 + 0.25 * Dw_equal);

fprintf('\nCalibration-only WMaxTwin scale weights:\n');
fprintf('  scale j     weight\n');
for ii = 1:numel(weights)
    fprintf('  %7d  %9.6f\n', weights(ii).j, weights(ii).weight);
end

%% Plot spectra and scale weights before the simulation
makeIntroFigures(X_cal_raw, wl128, y_cal, weights, OUTDIR);

%% Run repeated internal split comparison
% For each repetition, we construct the validation set inside the 172-sample
% calibration set.  The model is selected by validation MSE, then refit on all
% 172 calibration samples and evaluated once on the external test set.

rng(SEED);
methodNames = { ...
    'Random', 'SPlit', 'Twinning', 'DUPLEX', 'MaxTwin', ...
    'WMaxTwin gamma=0.10', 'WMaxTwin gamma=0.25', ...
    'WMaxTwin gamma=0.50', 'WMaxTwin gamma=0.75', ...
    'WMaxTwin equal gamma=0.25'};

results = struct('rep',{},'method',{},'jmax',{},'alpha',{},'val_mse',{},'test_mse',{},'test_rmse',{});
valCurves = struct('method',{},'jmax',{},'alpha',{},'val_mse',{});
resCount = 0;
vcCount = 0;

n_cal = size(X_cal,1);
allIdx = 1:n_cal;

for r = 0:(MATLAB_REPS - 1)
    rng(SEED + r);

    methods = struct();
    methods(1).name = 'Random';                    methods(1).val = randomSplit(n_cal, M_VAL);
    methods(2).name = 'SPlit';                     methods(2).val = supportSplit(D0, M_VAL);
    methods(3).name = 'Twinning';                  methods(3).val = twinningSplit(D0, M_VAL);
    methods(4).name = 'DUPLEX';                    methods(4).val = duplexSplit(D0, M_VAL);
    methods(5).name = 'MaxTwin';                   methods(5).val = maxtwinSplit(D0, M_VAL);
    methods(6).name = 'WMaxTwin gamma=0.10';       methods(6).val = maxtwinSplit(D_gamma{2}, M_VAL);
    methods(7).name = 'WMaxTwin gamma=0.25';       methods(7).val = maxtwinSplit(D_gamma{3}, M_VAL);
    methods(8).name = 'WMaxTwin gamma=0.50';       methods(8).val = maxtwinSplit(D_gamma{4}, M_VAL);
    methods(9).name = 'WMaxTwin gamma=0.75';       methods(9).val = maxtwinSplit(D_gamma{5}, M_VAL);
    methods(10).name = 'WMaxTwin equal gamma=0.25'; methods(10).val = maxtwinSplit(D_gamma_equal, M_VAL);

    for mm = 1:numel(methods)
        val_idx = unique(sort(methods(mm).val(:)'));
        if numel(val_idx) < M_VAL
            missing = setdiff(allIdx, val_idx);
            extra = missing(randperm(numel(missing), M_VAL - numel(val_idx)));
            val_idx = sort([val_idx, extra]);
        end
        train_idx = setdiff(allIdx, val_idx);

        res = fitSelectRefit(C_cal, y_cal, C_test, y_test, blocks, train_idx, val_idx, ALPHAS, JmaxGrid);

        resCount = resCount + 1;
        results(resCount).rep = r;
        results(resCount).method = methods(mm).name;
        results(resCount).jmax = res.jmax;
        results(resCount).alpha = res.alpha;
        results(resCount).val_mse = res.val_mse;
        results(resCount).test_mse = res.test_mse;
        results(resCount).test_rmse = res.test_rmse;

        if r == 0 && any(strcmp(methods(mm).name, {'Random','MaxTwin','WMaxTwin gamma=0.25','WMaxTwin gamma=0.50'}))
            curve = res.val_curve;
            for cc = 1:numel(curve)
                vcCount = vcCount + 1;
                valCurves(vcCount).method = methods(mm).name;
                valCurves(vcCount).jmax = curve(cc).jmax;
                valCurves(vcCount).alpha = curve(cc).alpha;
                valCurves(vcCount).val_mse = curve(cc).val_mse;
            end
        end
    end

    fprintf('Completed repetition %d / %d.\n', r+1, MATLAB_REPS);
end

%% Summarize, write tables, and draw figures
summary = summarizeResults(results, methodNames);
writeResultsTables(results, summary, weights, OUTDIR);
makeSimulationFigures(results, summary, valCurves, methodNames, OUTDIR);

fprintf('\nSummary table:\n');
disp(summary);
fprintf('\nOutputs written to:\n%s\n', OUTDIR);

%% Notes on interpretation
% MaxTwin corresponds to gamma = 0, that is, the raw standardized functional
% distance.  WMaxTwin uses gamma > 0 to introduce wavelet-scale information.
% This nested path is important because it shows that WMaxTwin is not a
% replacement for MaxTwin.  It is a controlled perturbation of MaxTwin toward a
% multiscale geometry.  In the paper, moderate gamma values improve external
% Tecator prediction relative to raw MaxTwin, while too much wavelet perturbation
% can again degrade performance.  This is the intended conclusion: WMaxTwin is
% useful when scale-local geometry is aligned with the validation problem, but it
% is not automatically superior for every weighting or every data set.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ensureTecatorData(path, url)
    if exist(path, 'file')
        info = dir(path);
        if info.bytes > 100000
            return;
        end
    end
    fprintf('Local Tecator file not found. Trying to download from StatLib...\n');
    try
        websave(path, url);
    catch ME
        error('Could not download Tecator data. Place tecator_statlib.txt in the working directory. Original error: %s', ME.message);
    end
end

function [X, endpoints, wavelengths] = loadTecator(path)
    txt = fileread(path);
    marker = 'extrapolation_examples=25';
    hits = strfind(txt, marker);
    if numel(hits) < 2
        error('Could not locate Tecator data block.');
    end
    start0 = hits(2);
    newlineRel = strfind(txt(start0:end), newline);
    if isempty(newlineRel)
        error('Could not locate numeric data after Tecator marker.');
    end
    startNumeric = start0 + newlineRel(1);
    nums = sscanf(txt(startNumeric:end), '%f');
    needed = 240 * 125;
    if numel(nums) < needed
        error('Parsed only %d numeric values; expected at least %d.', numel(nums), needed);
    end
    arr = reshape(nums(1:needed), 125, 240)';
    X = arr(:, 1:100);
    endpoints = arr(:, 123:125);  % moisture, fat, protein
    wavelengths = linspace(850, 1050, 100);
end

function [Xnew, grid] = interpToPower2(X, wavelengths, m)
    grid = linspace(min(wavelengths), max(wavelengths), m);
    Xnew = zeros(size(X,1), m);
    for i = 1:size(X,1)
        Xnew(i,:) = interp1(wavelengths, X(i,:), grid, 'linear', 'extrap');
    end
end

function sc = standardizeFit(X)
    sc.center = mean(X, 1);
    sc.scale = std(X, 0, 1);
    sc.scale(sc.scale < 1e-12) = 1;
end

function Xs = standardizeApply(X, sc)
    Xs = (X - sc.center) ./ sc.scale;
end

function [H, blocks] = haarTransformMatrix(n)
    J = log2(n);
    if abs(J - round(J)) > 1e-12
        error('n must be a power of two.');
    end
    J = round(J);
    H = zeros(n, n);
    row = 1;
    H(row,:) = ones(1,n) / sqrt(n);
    row = row + 1;
    blocks = struct();
    blocks(1).j = 0;
    blocks(1).idx = 1;
    bcount = 1;
    for lev = 1:J
        blockLen = n / 2^(lev - 1);
        half = blockLen / 2;
        nblocks = 2^(lev - 1);
        idxStart = row;
        for b = 1:nblocks
            v = zeros(1,n);
            start = (b - 1) * blockLen + 1;
            v(start:(start + half - 1)) = 1 / sqrt(blockLen);
            v((start + half):(start + blockLen - 1)) = -1 / sqrt(blockLen);
            H(row,:) = v;
            row = row + 1;
        end
        bcount = bcount + 1;
        blocks(bcount).j = lev;
        blocks(bcount).idx = idxStart:(row - 1);
    end
end

function D2 = pairwiseSqDists(Z)
    G = Z * Z';
    s = sum(Z.^2, 2);
    D2 = s + s' - 2 * G;
    D2(D2 < 0) = 0;
end

function D = normalizeDistance(D2)
    vals = D2(triu(true(size(D2)), 1));
    pos = vals(vals > 0);
    if isempty(pos)
        med = 1;
    else
        med = median(pos);
    end
    D = D2 / med;
end

function weights = scaleRelevanceWeights(C_cal, y_cal, blocks)
    y = (y_cal - mean(y_cal)) / std(y_cal, 0);
    vals = [];
    js = [];
    for bb = 1:numel(blocks)
        j = blocks(bb).j;
        if j == 0
            continue;
        end
        idx = blocks(bb).idx;
        Z = C_cal(:, idx);
        Zc = Z - mean(Z, 1);
        zs = std(Zc, 0, 1);
        ok = zs > 1e-12;
        if ~any(ok)
            score = 0;
        else
            Zok = Zc(:, ok);
            zsd = zs(ok);
            corr = (Zok' * y) ./ ((numel(y) - 1) * zsd');
            score = mean(corr.^2);
        end
        js(end+1) = j; %#ok<AGROW>
        vals(end+1) = score; %#ok<AGROW>
    end
    total = sum(vals);
    if total <= 0
        vals(:) = 1 / numel(vals);
    else
        vals = vals / total;
    end
    weights = struct('j', num2cell(js), 'weight', num2cell(vals));
end

function D = waveletWeightedDistance(C, blocks, weights)
    n = size(C,1);
    D = zeros(n,n);
    for ii = 1:numel(weights)
        j = weights(ii).j;
        idx = [];
        for bb = 1:numel(blocks)
            if blocks(bb).j == j
                idx = blocks(bb).idx;
                break;
            end
        end
        Z = C(:, idx);
        Zc = Z - mean(Z, 1);
        sdv = std(Zc, 0, 1);
        sdv(sdv < 1e-8) = 1;
        Zs = Zc ./ sdv;
        D = D + weights(ii).weight * normalizeDistance(pairwiseSqDists(Zs));
    end
end

function val = randomSplit(n, m_val)
    p = randperm(n);
    val = sort(p(1:m_val));
end

function val = supportSplit(D, m_val)
    n = size(D,1);
    Dj = D + 1e-10 * randn(n,n);
    meanD = mean(Dj, 2);
    selected = [];
    remaining = 1:n;
    sumMean = 0;
    sumPair = 0;
    for step = 0:(m_val - 1)
        bestObj = Inf;
        bestI = remaining(1);
        mNew = step + 1;
        for rr = 1:numel(remaining)
            i = remaining(rr);
            if step == 0
                addPair = 0;
            else
                addPair = 2 * sum(Dj(i, selected));
            end
            obj = 2 * (sumMean + meanD(i)) / mNew - (sumPair + addPair) / (mNew^2);
            if obj < bestObj
                bestObj = obj;
                bestI = i;
            end
        end
        if step > 0
            sumPair = sumPair + 2 * sum(Dj(bestI, selected));
        end
        selected = [selected, bestI]; %#ok<AGROW>
        remaining(remaining == bestI) = [];
        sumMean = sumMean + meanD(bestI);
    end
    val = selected;
end

function val = duplexSplit(D, m_val)
    n = size(D,1);
    Dj = D + 1e-10 * randn(n,n);
    U = triu(Dj, 1);
    [~, linear] = max(U(:));
    [i,j] = ind2sub(size(U), linear);
    selected = [i, j];
    remaining = setdiff(1:n, selected);
    while numel(selected) < m_val
        mind = zeros(numel(remaining),1);
        for rr = 1:numel(remaining)
            mind(rr) = min(Dj(remaining(rr), selected));
        end
        [~,loc] = max(mind);
        chosen = remaining(loc);
        selected = [selected, chosen]; %#ok<AGROW>
        remaining(remaining == chosen) = [];
    end
    val = selected;
end

function val = twinningSplit(D, m_val)
    n = size(D,1);
    anchors = supportSplit(D, m_val);
    Dj = D + 1e-10 * randn(n,n);
    unused = 1:n;
    val = [];
    for aa = 1:numel(anchors)
        if numel(val) >= m_val
            break;
        end
        a = anchors(aa);
        if ~ismember(a, unused) || numel(unused) < 2
            continue;
        end
        candidates = setdiff(unused, a);
        [~,loc] = min(Dj(a, candidates));
        b = candidates(loc);
        if rand < 0.5
            val(end+1) = a; %#ok<AGROW>
        else
            val(end+1) = b; %#ok<AGROW>
        end
        unused = setdiff(unused, [a, b]);
    end
    while numel(val) < m_val && numel(unused) >= 2
        subD = Dj(unused, unused);
        subD(1:size(subD,1)+1:end) = Inf;
        [~,linear] = min(subD(:));
        [ia, ib] = ind2sub(size(subD), linear);
        a = unused(ia); b = unused(ib);
        if rand < 0.5
            val(end+1) = a; %#ok<AGROW>
        else
            val(end+1) = b; %#ok<AGROW>
        end
        unused = setdiff(unused, [a, b]);
    end
    if numel(val) < m_val
        p = randperm(numel(unused));
        val = [val, unused(p(1:(m_val - numel(val))))];
    end
end

function val = maxtwinSplit(D, m_val)
    n = size(D,1);
    Dj = D + 1e-10 * randn(n,n);
    unused = 1:n;
    val = [];
    while numel(val) < m_val && numel(unused) >= 2
        subD = Dj(unused, unused);
        subD(1:size(subD,1)+1:end) = Inf;
        [~,linear] = min(subD(:));
        [ia, ib] = ind2sub(size(subD), linear);
        a = unused(ia); b = unused(ib);
        if rand < 0.5
            val(end+1) = a; %#ok<AGROW>
        else
            val(end+1) = b; %#ok<AGROW>
        end
        unused = setdiff(unused, [a, b]);
    end
    if numel(val) < m_val
        p = randperm(numel(unused));
        val = [val, unused(p(1:(m_val - numel(val))))];
    end
end

function F = featuresByJmax(C, blocks, jmax)
    idx = 1;
    for bb = 1:numel(blocks)
        if blocks(bb).j >= 1 && blocks(bb).j <= jmax
            idx = [idx, blocks(bb).idx]; %#ok<AGROW>
        end
    end
    F = C(:, idx);
end

function fit = ridgeFit(X, y, alpha)
    Xaug = [ones(size(X,1),1), X];
    p = size(Xaug,2);
    P = diag([0; alpha * ones(p - 1,1)]);
    beta = (Xaug' * Xaug + P) \ (Xaug' * y);
    fit.beta = beta;
end

function yhat = ridgePredict(fit, X)
    Xaug = [ones(size(X,1),1), X];
    yhat = Xaug * fit.beta;
end

function res = fitSelectRefit(C_cal, y_cal, C_test, y_test, blocks, train_idx, val_idx, alphas, jmax_grid)
    bestVmse = Inf;
    bestJ = NaN;
    bestAlpha = NaN;
    curve = struct('jmax',{},'alpha',{},'val_mse',{});
    cc = 0;
    for jj = 1:numel(jmax_grid)
        jmax = jmax_grid(jj);
        F = featuresByJmax(C_cal, blocks, jmax);
        Xtr = F(train_idx, :);
        Xva = F(val_idx, :);
        ytr = y_cal(train_idx);
        yva = y_cal(val_idx);
        sc = standardizeFit(Xtr);
        Xtr_s = standardizeApply(Xtr, sc);
        Xva_s = standardizeApply(Xva, sc);
        for aa = 1:numel(alphas)
            alpha = alphas(aa);
            fit = ridgeFit(Xtr_s, ytr, alpha);
            pred = ridgePredict(fit, Xva_s);
            vmse = mean((yva - pred).^2);
            cc = cc + 1;
            curve(cc).jmax = jmax;
            curve(cc).alpha = alpha;
            curve(cc).val_mse = vmse;
            if vmse < bestVmse
                bestVmse = vmse;
                bestJ = jmax;
                bestAlpha = alpha;
            end
        end
    end
    Fcal = featuresByJmax(C_cal, blocks, bestJ);
    Ftest = featuresByJmax(C_test, blocks, bestJ);
    sc = standardizeFit(Fcal);
    fit = ridgeFit(standardizeApply(Fcal, sc), y_cal, bestAlpha);
    predTest = ridgePredict(fit, standardizeApply(Ftest, sc));
    testMSE = mean((y_test - predTest).^2);
    res.jmax = bestJ;
    res.alpha = bestAlpha;
    res.val_mse = bestVmse;
    res.test_mse = testMSE;
    res.test_rmse = sqrt(testMSE);
    res.pred_test = predTest;
    res.val_curve = curve;
end

function summary = summarizeResults(results, methodNames)
    rows = {};
    for mm = 1:numel(methodNames)
        name = methodNames{mm};
        idx = strcmp({results.method}, name);
        if ~any(idx), continue; end
        sub = results(idx);
        jmax = [sub.jmax]';
        alpha = [sub.alpha]';
        valm = [sub.val_mse]';
        rmse = [sub.test_rmse]';
        mse = [sub.test_mse]';
        rows(end+1,:) = {name, mean(jmax), localStd(jmax), median(alpha), mean(valm), localStd(valm), mean(rmse), localStd(rmse), mean(mse), localStd(mse)}; %#ok<AGROW>
    end
    summary = cell2table(rows, 'VariableNames', {'method','mean_jmax','sd_jmax','median_alpha','mean_val_mse','sd_val_mse','mean_test_rmse','sd_test_rmse','mean_test_mse','sd_test_mse'});
end

function s = localStd(x)
    if numel(x) <= 1
        s = 0;
    else
        s = std(x, 0);
    end
end

function writeResultsTables(results, summary, weights, outdir)
    n = numel(results);
    rep = [results.rep]';
    method = {results.method}';
    jmax = [results.jmax]';
    alpha = [results.alpha]';
    val_mse = [results.val_mse]';
    test_mse = [results.test_mse]';
    test_rmse = [results.test_rmse]';
    T = table(rep, method, jmax, alpha, val_mse, test_mse, test_rmse);
    writetable(T, fullfile(outdir, 'tecator_split_results_MATLAB.csv'));
    writetable(summary, fullfile(outdir, 'tecator_split_summary_MATLAB.csv'));
    scale_j = [weights.j]';
    response_weight = [weights.weight]';
    equal_weight = ones(size(response_weight)) / numel(response_weight);
    WT = table(scale_j, response_weight, equal_weight);
    writetable(WT, fullfile(outdir, 'tecator_scale_weights_MATLAB.csv'));

    fid = fopen(fullfile(outdir, 'tecator_wmaxtwin_table_MATLAB.tex'), 'w');
    fprintf(fid, '%% Generated by WMaxTwin_Tecator_Functional_Regression_MATLAB_Notebook.m\n\n');
    fprintf(fid, '\\begin{tabular}{lrrrrr}\n');
    fprintf(fid, '\\hline\n');
    fprintf(fid, 'split & mean $\\widehat{j}_{\\max}$ & median $\\widehat\\lambda$ & mean val. MSE & mean test RMSE & sd test RMSE \\\\\n');
    fprintf(fid, '\\hline\n');
    for i = 1:height(summary)
        fprintf(fid, '%s & %.2f & %.3g & %.4f & %.4f & %.4f \\\\\n', ...
            summary.method{i}, summary.mean_jmax(i), summary.median_alpha(i), summary.mean_val_mse(i), summary.mean_test_rmse(i), summary.sd_test_rmse(i));
    end
    fprintf(fid, '\\hline\n');
    fprintf(fid, '\\end{tabular}\n');
    fclose(fid);
end

function makeIntroFigures(X_cal_raw, wl, y_cal, weights, outdir)
    f = figure('Visible','off', 'Position', [100 100 900 520]);
    hold on;
    for i = 1:4:size(X_cal_raw,1)
        plot(wl, X_cal_raw(i,:), 'Color', [0.70 0.70 0.70]);
    end
    xlabel('Wavelength (nm)'); ylabel('Absorbance');
    title('Tecator absorbance spectra, calibration samples');
    box on; hold off;
    print(f, fullfile(outdir, 'fig1_tecator_spectra.png'), '-dpng', '-r220');
    close(f);

    [~,ord] = sort(y_cal);
    n = numel(y_cal);
    group = zeros(n,1);
    for q = 1:4
        lo = floor((q-1)*n/4) + 1;
        hi = floor(q*n/4);
        group(ord(lo:hi)) = q;
    end
    f = figure('Visible','off', 'Position', [100 100 900 520]);
    hold on;
    styles = {'-','--',':','-.'};
    for q = 1:4
        m = mean(X_cal_raw(group == q, :), 1);
        plot(wl, m, 'LineWidth', 2, 'LineStyle', styles{q});
    end
    xlabel('Wavelength (nm)'); ylabel('Mean absorbance');
    title('Mean spectra by fat-content quartile');
    legend({'fat quartile 1','fat quartile 2','fat quartile 3','fat quartile 4'}, 'Location','best', 'Box','off');
    box on; hold off;
    print(f, fullfile(outdir, 'fig2_quartile_mean_spectra.png'), '-dpng', '-r220');
    close(f);

    f = figure('Visible','off', 'Position', [100 100 760 500]);
    bar([weights.j], [weights.weight]);
    xlabel('Wavelet scale j (coarse to fine)');
    ylabel('Normalized response-relevance weight');
    title('Calibration-only scale weights for WMaxTwin');
    box on;
    print(f, fullfile(outdir, 'fig3_response_scale_weights.png'), '-dpng', '-r220');
    close(f);
end

function makeSimulationFigures(results, summary, valCurves, methodNames, outdir)
    methods = summary.method;

    if ~isempty(valCurves)
        f = figure('Visible','off', 'Position', [100 100 950 560]);
        hold on;
        labs = unique({valCurves.method}, 'stable');
        styles = {'-o','-s','-^','-d','-*'};
        for ii = 1:numel(labs)
            idx = strcmp({valCurves.method}, labs{ii});
            sub = valCurves(idx);
            jvals = unique([sub.jmax]);
            y = zeros(size(jvals));
            for jj = 1:numel(jvals)
                y(jj) = min([sub([sub.jmax] == jvals(jj)).val_mse]);
            end
            plot(jvals, y, styles{1 + mod(ii-1,numel(styles))}, 'LineWidth', 1.5, 'MarkerSize', 5);
        end
        xlabel('Maximum included wavelet scale jmax');
        ylabel('Best validation MSE over ridge penalties');
        title('Example validation curves, first repetition');
        legend(labs, 'Location','best', 'Box','off');
        box on; hold off;
        print(f, fullfile(outdir, 'fig4_validation_curves_rep0.png'), '-dpng', '-r220');
        close(f);
    end

    f = figure('Visible','off', 'Position', [100 100 1200 650]);
    vals = cell(numel(methods),1);
    for i = 1:numel(methods)
        vals{i} = [results(strcmp({results.method}, methods{i})).test_rmse];
    end
    drawSimpleBoxplot(vals, methods, 'External test RMSE: fat percent', 'Tecator external test error across repeated internal splits');
    print(f, fullfile(outdir, 'fig5_test_rmse_boxplot.png'), '-dpng', '-r220');
    close(f);

    f = figure('Visible','off', 'Position', [100 100 1200 650]);
    vals = cell(numel(methods),1);
    for i = 1:numel(methods)
        vals{i} = [results(strcmp({results.method}, methods{i})).jmax];
    end
    drawSimpleBoxplot(vals, methods, 'Selected jmax', 'Selected wavelet-resolution cutoff');
    print(f, fullfile(outdir, 'fig6_selected_jmax_boxplot.png'), '-dpng', '-r220');
    close(f);

    f = figure('Visible','off', 'Position', [100 100 1200 650]);
    vals = cell(numel(methods),1);
    for i = 1:numel(methods)
        vals{i} = log10([results(strcmp({results.method}, methods{i})).alpha]);
    end
    drawSimpleBoxplot(vals, methods, 'log10 selected ridge penalty', 'Selected shrinkage penalty across repeated internal splits');
    print(f, fullfile(outdir, 'fig7_selected_alpha_boxplot.png'), '-dpng', '-r220');
    close(f);

    gammaLabels = {'MaxTwin','WMaxTwin gamma=0.10','WMaxTwin gamma=0.25','WMaxTwin gamma=0.50','WMaxTwin gamma=0.75'};
    xvals = [0, 0.10, 0.25, 0.50, 0.75];
    y = nan(size(xvals)); e = nan(size(xvals));
    for i = 1:numel(gammaLabels)
        idx = strcmp({results.method}, gammaLabels{i});
        v = [results(idx).test_rmse];
        if ~isempty(v)
            y(i) = mean(v);
            if numel(v) > 1
                e(i) = std(v,0) / sqrt(numel(v));
            else
                e(i) = 0;
            end
        end
    end
    f = figure('Visible','off', 'Position', [100 100 850 560]);
    errorbar(xvals, y, e, '-o', 'LineWidth', 1.8, 'MarkerSize', 6);
    xlabel('Nested WMaxTwin mixing parameter gamma');
    ylabel('Mean external test RMSE');
    title('Nested path: MaxTwin is gamma=0');
    box on;
    print(f, fullfile(outdir, 'fig8_gamma_path_rmse.png'), '-dpng', '-r220');
    close(f);
end

function drawSimpleBoxplot(vals, labels, ylab, ttl)
    hold on;
    n = numel(vals);
    for i = 1:n
        x = vals{i};
        x = x(~isnan(x));
        q1 = localPercentile(x, 25);
        q2 = localPercentile(x, 50);
        q3 = localPercentile(x, 75);
        iqr = q3 - q1;
        lo = max(min(x), q1 - 1.5 * iqr);
        hi = min(max(x), q3 + 1.5 * iqr);
        rectangle('Position', [i - 0.25, q1, 0.5, max(q3 - q1, eps)], 'EdgeColor', 'k');
        line([i - 0.25, i + 0.25], [q2, q2], 'Color','k', 'LineWidth', 1.5);
        line([i, i], [lo, q1], 'Color','k');
        line([i, i], [q3, hi], 'Color','k');
        line([i - 0.15, i + 0.15], [lo, lo], 'Color','k');
        line([i - 0.15, i + 0.15], [hi, hi], 'Color','k');
        plot(i, mean(x), 'ko', 'MarkerFaceColor','k', 'MarkerSize', 4);
    end
    xlim([0.5, n + 0.5]);
    set(gca, 'XTick', 1:n, 'XTickLabel', labels, 'XTickLabelRotation', 35);
    ylabel(ylab); title(ttl); box on; hold off;
end

function q = localPercentile(x, p)
    x = sort(x(:));
    if isempty(x), q = NaN; return; end
    if numel(x) == 1, q = x; return; end
    pos = 1 + (p/100) * (numel(x) - 1);
    lo = floor(pos); hi = ceil(pos);
    if lo == hi
        q = x(lo);
    else
        q = x(lo) + (pos - lo) * (x(hi) - x(lo));
    end
end
