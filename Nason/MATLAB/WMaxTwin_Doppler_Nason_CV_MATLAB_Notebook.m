%% WMaxTwin Doppler / Nason-type wavelet cross-validation in MATLAB
% Multiscale MaxTwin: Example II
%
% This is a self-contained MATLAB translation of the Python and R notebooks
% for the Doppler / Nason-type wavelet cross-validation experiment.
%
% Main statistical point:
%   Nason's even/odd split is grid-compatible.
%   MaxTwin and WMaxTwin masks are generally not dyadic subgrids.
%   Therefore the fitting step is written as a masked wavelet penalized
%   least-squares problem on the original grid, not as an ordinary DWT of
%   the retained observations.
%
% Dependencies:
%   Base MATLAB only.  The periodic DWT/IDWT for sym4 and db2 is implemented
%   directly from filter coefficients.  No Wavelet Toolbox is required.
%
% GitHub note:
%   This sectioned .m file is intended to behave like a MATLAB notebook.
%   Open it in MATLAB and run section-by-section using the %% cells.  It can
%   also be saved as a MATLAB Live Script (.mlx) from the MATLAB Editor.

clear; close all; clc;

%% User settings
N          = 1024;      % signal length, power of two
SNR        = 7;         % signal-to-noise ratio convention used in the paper
SIGMA      = 1;         % noise standard deviation
WAVELET    = 'sym4';    % 'sym4' = 8-tap Symmlet; 'db2' also provided
J0         = 3;         % coarsest detail level kept in the thresholding family
NREP       = 18;        % default reproduces the paper/notebook comparison
SEED       = 2026;
N_ITER     = 30;        % FISTA iterations for each masked fit
WEIGHT_MODE = 'paper';  % 'paper' uses explicit omega_j; 'script' reproduces original code weights
OUTPUT_DIR = fullfile(pwd, 'outputs_matlab');

if ~exist(OUTPUT_DIR, 'dir')
    mkdir(OUTPUT_DIR);
end

rng(SEED, 'twister');

%% Statistical setup
% We generate the Doppler signal on the original dyadic grid and add Gaussian
% noise.  The target is the true signal.  The competing procedures differ in
% how they select thresholds and how they construct the training/validation
% split for cross-validation.
%
% Methods compared:
%   1. Noisy data.
%   2. Levelwise SURE benchmark.
%   3. Nason even/odd CV with a fixed global threshold.
%   4. Location-only MaxTwin with a fixed global threshold.
%   5. WMaxTwin with high-frequency feature geometry and a ramp threshold.
%   6. Oracle ramp, included only as an unattainable benchmark.
%
% WMaxTwin does not itself produce thresholds.  It constructs a scale-aware
% validation geometry.  The thresholding rule is a separate estimator-side
% choice.  In this example, WMaxTwin is paired with a one-parameter ramp
% threshold family because the Doppler signal has difficult high-frequency
% oscillatory behavior.

[t0, f0] = doppler_signal(N, SNR, SIGMA);
J = log2(N);
weights = wmaxtwin_scale_weights(J0, J, WEIGHT_MODE);

disp('WMaxTwin feature weights for levels j = J0:(J-1):');
disp(table((J0:J-1)', weights(:), 'VariableNames', {'Level_j','omega_j'}));

%% Figure 1: clean Doppler signal and one noisy realization
rng(SEED, 'twister');
y0 = f0 + SIGMA * randn(N, 1);

fig = figure('Name', 'Doppler signal');
plot(t0, y0, 'Color', [0.75 0.75 0.75], 'LineWidth', 0.7); hold on;
plot(t0, f0, 'k', 'LineWidth', 2.0);
xlabel('t'); ylabel('signal');
title('Doppler signal: one noisy realization');
legend({'noisy observation','true signal'}, 'Location', 'best');
save_figure(fig, fullfile(OUTPUT_DIR, 'fig01_doppler_signal.png'));

%% Figure 2: WMaxTwin feature weights
fig = figure('Name', 'WMaxTwin feature weights');
bar(J0:J-1, weights);
xlabel('wavelet level j'); ylabel('\omega_j');
title('WMaxTwin feature weights used in the split geometry');
save_figure(fig, fullfile(OUTPUT_DIR, 'fig02_feature_weights.png'));

%% Run the Monte Carlo simulation
% The default NREP = 18 gives the same scale of calculation as the Python
% notebook.  For a fast smoke test set NREP = 2 or NREP = 3.  For a more stable
% paper-scale run, increase NREP.

result = run_simulation(N, SNR, SIGMA, WAVELET, J0, NREP, SEED, N_ITER, WEIGHT_MODE, OUTPUT_DIR);
summary_table = result.summary;
records = result.records;
saved = result.saved;

disp(summary_table);

%% Interpretation of the Monte Carlo table
% In the intended run, the WMaxTwin ramp method should be close to the oracle
% ramp and below the fixed-threshold MaxTwin method.  This is the substantive
% point of the example.  It does not say that WMaxTwin alone is a universal
% improvement.  The improvement comes from the coherent combination
%
%   WMaxTwin split geometry + masked wavelet fit + high-frequency ramp family.
%
% The MaxTwin fixed-threshold method changes the split but keeps a less suitable
% threshold family.  It is therefore not expected to dominate the Nason benchmark.

%% Local functions
% MATLAB permits local functions at the end of a script.  This keeps the file
% self-contained and GitHub-friendly.

function [h, g] = wavelet_filters(wavelet)
    switch lower(wavelet)
        case 'sym4'
            h = [
                -0.07576571478927333;
                -0.02963552764599851;
                 0.49761866763201545;
                 0.8037387518059161;
                 0.29785779560527736;
                -0.09921954357684722;
                -0.012603967262037833;
                 0.0322231006040427
            ];
            g = [
                -0.0322231006040427;
                -0.012603967262037833;
                 0.09921954357684722;
                 0.29785779560527736;
                -0.8037387518059161;
                 0.49761866763201545;
                 0.02963552764599851;
                -0.07576571478927333
            ];
        case 'db2'
            h = [
                -0.12940952255126034;
                 0.2241438680420134;
                 0.8365163037378079;
                 0.48296291314453416
            ];
            g = [
                -0.48296291314453416;
                 0.8365163037378079;
                -0.2241438680420134;
                -0.12940952255126034
            ];
        otherwise
            error('Unknown wavelet. Use ''sym4'' or ''db2''.');
    end
end

function idx = periodic_idx(n, L)
    m = n / 2;
    row_base = 2 * (0:(m - 1))';
    offsets = 0:(L - 1);
    idx = mod(row_base + offsets, n) + 1;
end

function [a, d] = dwt_step(x, h, g)
    x = x(:);
    idx = periodic_idx(length(x), length(h));
    X = x(idx);
    a = X * h(:);
    d = X * g(:);
end

function x = idwt_step(a, d, h, g)
    a = a(:);
    d = d(:);
    m = length(a);
    n = 2 * m;
    idx = periodic_idx(n, length(h));
    vals = a * h(:)' + d * g(:)';
    x = accumarray(idx(:), vals(:), [n, 1], @sum, 0);
end

function dec = wavedec_periodic(x, j0, wavelet)
    x = x(:);
    [h, g] = wavelet_filters(wavelet);
    n = length(x);
    J = log2(n);
    if abs(J - round(J)) > 1e-12
        error('n must be a power of two.');
    end
    J = round(J);
    if ~(1 <= j0 && j0 < J)
        error('j0 must satisfy 1 <= j0 < J.');
    end

    a = x;
    details_by_level = cell(J, 1);
    for lev = (J - 1):-1:j0
        [a, d] = dwt_step(a, h, g);
        details_by_level{lev + 1} = d;
    end

    details = cell(J - j0, 1);
    for r = 1:(J - j0)
        lev = j0 + r - 1;
        details{r} = details_by_level{lev + 1};
    end

    dec.a = a;
    dec.details = details;
    dec.levels = (j0:(J - 1))';
end

function x = waverec_periodic(a, details, j0, wavelet)
    %#ok<INUSD> j0 kept for readability and consistency with wavedec_periodic.
    [h, g] = wavelet_filters(wavelet);
    cur = a(:);
    for r = 1:length(details)
        cur = idwt_step(cur, details{r}, h, g);
    end
    x = cur;
end

function y = soft_thresh(x, lambda)
    y = sign(x) .* max(abs(x) - lambda, 0);
end

function fit = shrink_full(y, lams, j0, wavelet)
    dec = wavedec_periodic(y, j0, wavelet);
    details2 = cell(size(dec.details));
    for r = 1:length(dec.details)
        details2{r} = soft_thresh(dec.details{r}, lams(r));
    end
    fit = waverec_periodic(dec.a, details2, j0, wavelet);
end

function [t, f] = doppler_signal(n, snr, sigma)
    t = ((1:n)' - 0.5) / n;
    f = sqrt(t .* (1 - t)) .* sin(2 * pi * 1.05 ./ (t + 0.05));
    f = f - mean(f);
    f = f * sqrt(snr * sigma^2 / var(f, 1));
end

function val = mse(a, b)
    val = mean((a(:) - b(:)).^2);
end

function lambda = sure_threshold_detail(d, sigma)
    x = sort(abs(d(:)));
    n = length(x);
    s2 = sigma^2;
    csum = cumsum(x.^2);
    k = (1:n)';
    risks = n * s2 + csum + (n - k) .* x.^2 - 2 * s2 * k;
    [~, pos] = min(risks);
    lambda = x(pos);
end

function lams = level_sure_lams(y, j0, wavelet, sigma)
    J = round(log2(length(y)));
    dec = wavedec_periodic(y, j0, wavelet);
    lams = zeros(J - j0, 1);
    for r = 1:length(dec.details)
        lams(r) = sure_threshold_detail(dec.details{r}, sigma);
    end
end

function lams = fixed_lams(lambda, j0, J)
    lams = lambda * ones(J - j0, 1);
end

function lams = ramp_lams(slope, j0, J, power)
    js = (j0:(J - 1))';
    lams = slope * ((js - j0) / (J - 1 - j0)).^power;
end

function fit = fit_masked_wavelet(y, mask, lams, j0, wavelet, n_iter)
    y = y(:);
    mask = mask(:);
    dec = wavedec_periodic(y .* mask, j0, wavelet);
    a = dec.a;
    details = dec.details;
    za = a;
    zd = details;
    tk = 1;

    for iter = 1:n_iter %#ok<NASGU>
        fz = waverec_periodic(za, zd, j0, wavelet);
        resid = mask .* (fz - y);
        grad = wavedec_periodic(resid, j0, wavelet);

        a_new = za - grad.a;
        d_new = cell(size(zd));
        for r = 1:length(zd)
            d_new{r} = soft_thresh(zd{r} - grad.details{r}, lams(r));
        end

        t_new = (1 + sqrt(1 + 4 * tk^2)) / 2;
        beta = (tk - 1) / t_new;
        za = a_new + beta * (a_new - a);
        zd_new = cell(size(zd));
        for r = 1:length(zd)
            zd_new{r} = d_new{r} + beta * (d_new{r} - details{r});
        end

        a = a_new;
        details = d_new;
        zd = zd_new;
        tk = t_new;
    end

    fit = waverec_periodic(a, details, j0, wavelet);
end

function cont = detail_contributions(y, j0, wavelet)
    dec = wavedec_periodic(y, j0, wavelet);
    zero_a = zeros(size(dec.a));
    cont = cell(size(dec.details));
    for r = 1:length(dec.details)
        tmp = cell(size(dec.details));
        for s = 1:length(dec.details)
            tmp{s} = zeros(size(dec.details{s}));
        end
        tmp{r} = dec.details{r};
        cont{r} = waverec_periodic(zero_a, tmp, j0, wavelet);
    end
end

function weights = wmaxtwin_scale_weights(j0, J, mode)
    js = (j0:(J - 1))';
    if strcmpi(mode, 'paper') && j0 == 3 && J == 10
        weights = [0.015; 0.036; 0.074; 0.124; 0.182; 0.248; 0.321];
        return;
    end

    if strcmpi(mode, 'script')
        weights = zeros(length(js), 1);
        for ii = 1:length(js)
            j = js(ii);
            if j >= J - 3
                rel = j - (J - 3) + 1;
                tmp = [1.0, 1.8, 2.6];
                weights(ii) = tmp(rel);
            else
                weights(ii) = 0.10 * ((j - j0 + 1) / (J - j0))^2;
            end
        end
        return;
    end

    eta = 0.05;
    weights = eta + ((js - j0) / (J - 1 - j0)).^1.5;
    weights = weights / sum(weights);
end

function z = standardize_vec(x)
    x = x(:);
    sx = std(x, 1);
    if sx < 1e-12
        z = zeros(size(x));
    else
        z = (x - mean(x)) / sx;
    end
end

function split = max_or_wmax_split(y, j0, wavelet, kind, weight_mode)
    y = y(:);
    n = length(y);
    t = ((1:n)' - 0.5) / n;

    if strcmpi(kind, 'nason')
        split.train = (1:2:n)';
        split.val = (2:2:n)';
        return;
    end

    if strcmpi(kind, 'max')
        m = n / 2;
        train = zeros(m, 1);
        val = zeros(m, 1);
        for q = 1:m
            i = 2 * q - 1;
            j = 2 * q;
            if mod(q, 2) == 1
                train(q) = i; val(q) = j;
            else
                train(q) = j; val(q) = i;
            end
        end
        split.train = train;
        split.val = val;
        return;
    end

    cont = detail_contributions(y, j0, wavelet);
    J = round(log2(n));
    weights = wmaxtwin_scale_weights(j0, J, weight_mode);

    X = zeros(n, 1 + (J - j0));
    X(:, 1) = standardize_vec(t);
    for j = j0:(J - 1)
        r = j - j0 + 1;
        q = abs(cont{r});
        q = standardize_vec(q);
        X(:, r + 1) = sqrt(weights(r)) * q;
    end

    activity = zeros(n, 1);
    for j = max(j0, J - 3):(J - 1)
        r = j - j0 + 1;
        activity = activity + abs(cont{r});
    end
    [~, ord] = sort(activity, 'descend');

    remaining = true(n, 1);
    pairs = zeros(n / 2, 2);
    cnt = 0;
    for ii = 1:length(ord)
        i = ord(ii);
        if ~remaining(i)
            continue;
        end
        remaining(i) = false;
        rem = find(remaining);
        if isempty(rem)
            break;
        end
        diffs = X(rem, :) - X(i, :);
        dist = sum(diffs.^2, 2);
        [~, pos] = min(dist);
        jj = rem(pos);
        remaining(jj) = false;
        cnt = cnt + 1;
        pairs(cnt, :) = [i, jj];
    end
    pairs = pairs(1:cnt, :);

    bal = zeros(1, size(X, 2));
    train = zeros(cnt, 1);
    val = zeros(cnt, 1);
    for q = 1:cnt
        i = pairs(q, 1);
        j = pairs(q, 2);
        diff = X(i, :) - X(j, :);
        if norm(bal + diff) <= norm(bal - diff)
            train(q) = i; val(q) = j;
            bal = bal + diff;
        else
            train(q) = j; val(q) = i;
            bal = bal - diff;
        end
    end

    split.train = train;
    split.val = val;
end

function cv = cv_select_fixed(y, train, val, j0, wavelet, n_iter)
    n = length(y);
    J = round(log2(n));
    mask_t = zeros(n, 1); mask_t(train) = 1;
    mask_v = zeros(n, 1); mask_v(val) = 1;

    lambdas = linspace(0.2, 2.4, 12)';
    scores = zeros(length(lambdas), 1);
    for ii = 1:length(lambdas)
        lams = fixed_lams(lambdas(ii), j0, J);
        pred_t = fit_masked_wavelet(y, mask_t, lams, j0, wavelet, n_iter);
        pred_v = fit_masked_wavelet(y, mask_v, lams, j0, wavelet, n_iter);
        scores(ii) = mean((pred_t(val) - y(val)).^2) + mean((pred_v(train) - y(train)).^2);
    end

    rows = table(scores, lambdas, 'VariableNames', {'score','lambda'});
    rows = sortrows(rows, 'score');
    cv.rows = rows;
    cv.chosen = rows(1, :);
    cv.lams = fixed_lams(rows.lambda(1), j0, J);
end

function cv = cv_select_ramp(y, train, val, j0, wavelet, n_iter, flat_tol)
    n = length(y);
    J = round(log2(n));
    mask_t = zeros(n, 1); mask_t(train) = 1;
    mask_v = zeros(n, 1); mask_v(val) = 1;

    slopes = linspace(0.6, 3.4, 15)';
    scores = zeros(length(slopes), 1);
    for ii = 1:length(slopes)
        lams = ramp_lams(slopes(ii), j0, J, 1.5);
        pred_t = fit_masked_wavelet(y, mask_t, lams, j0, wavelet, n_iter);
        pred_v = fit_masked_wavelet(y, mask_v, lams, j0, wavelet, n_iter);
        scores(ii) = mean((pred_t(val) - y(val)).^2) + mean((pred_v(train) - y(train)).^2);
    end

    rows = table(scores, slopes, 'VariableNames', {'score','slope'});
    rows = sortrows(rows, 'score');
    min_score = rows.score(1);
    cand = rows(rows.score <= flat_tol * min_score, :);
    [~, pos] = max(cand.slope);
    cv.rows = rows;
    cv.chosen = cand(pos, :);
    cv.lams = ramp_lams(cv.chosen.slope(1), j0, J, 1.5);
end

function out = oracle_ramp(y, f, j0, wavelet)
    J = round(log2(length(y)));
    slopes = linspace(0.4, 4.0, 37);
    best_err = inf;
    best_slope = NaN;
    for ii = 1:length(slopes)
        lams = ramp_lams(slopes(ii), j0, J, 1.5);
        err = mse(shrink_full(y, lams, j0, wavelet), f);
        if err < best_err
            best_err = err;
            best_slope = slopes(ii);
        end
    end
    out.error = best_err;
    out.slope = best_slope;
end

function result = run_simulation(n, snr, sigma, wavelet, j0, nrep, seed, n_iter, weight_mode, output_dir)
    J = round(log2(n));
    if 2^J ~= n
        error('n must be a power of two.');
    end
    if j0 ~= 3
        warning('This example was tuned for j0 = 3.');
    end

    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end

    [t, f] = doppler_signal(n, snr, sigma);
    rng(seed, 'twister');

    rep_id = zeros(nrep, 1);
    Noisy = zeros(nrep, 1);
    LevelSure = zeros(nrep, 1);
    NasonEvenOddCV = zeros(nrep, 1);
    MaxTwinFixedCV = zeros(nrep, 1);
    WMaxTwinRampCV = zeros(nrep, 1);
    OracleRamp = zeros(nrep, 1);
    WMaxTwinSlope = zeros(nrep, 1);
    OracleSlope = zeros(nrep, 1);

    saved = struct();
    saved_exists = false;

    for r = 1:nrep
        fprintf('replication %d of %d\n', r, nrep);
        y = f + sigma * randn(n, 1);

        l_sure = level_sure_lams(y, j0, wavelet, sigma);
        f_sure = shrink_full(y, l_sure, j0, wavelet);

        sp_n = max_or_wmax_split(y, j0, wavelet, 'nason', weight_mode);
        cv_n = cv_select_fixed(y, sp_n.train, sp_n.val, j0, wavelet, n_iter);
        f_n = shrink_full(y, cv_n.lams, j0, wavelet);

        sp_m = max_or_wmax_split(y, j0, wavelet, 'max', weight_mode);
        cv_m = cv_select_fixed(y, sp_m.train, sp_m.val, j0, wavelet, n_iter);
        f_m = shrink_full(y, cv_m.lams, j0, wavelet);

        sp_w = max_or_wmax_split(y, j0, wavelet, 'wmax', weight_mode);
        cv_w = cv_select_ramp(y, sp_w.train, sp_w.val, j0, wavelet, n_iter, 1.015);
        f_w = shrink_full(y, cv_w.lams, j0, wavelet);

        or = oracle_ramp(y, f, j0, wavelet);

        rep_id(r) = r;
        Noisy(r) = mse(y, f);
        LevelSure(r) = mse(f_sure, f);
        NasonEvenOddCV(r) = mse(f_n, f);
        MaxTwinFixedCV(r) = mse(f_m, f);
        WMaxTwinRampCV(r) = mse(f_w, f);
        OracleRamp(r) = or.error;
        WMaxTwinSlope(r) = cv_w.chosen.slope(1);
        OracleSlope(r) = or.slope;

        if ~saved_exists
            saved.t = t;
            saved.f = f;
            saved.y = y;
            saved.f_sure = f_sure;
            saved.f_nason = f_n;
            saved.f_maxtwin = f_m;
            saved.f_wmax = f_w;
            saved.l_sure = l_sure;
            saved.l_wmax = cv_w.lams;
            saved.rows_w = cv_w.rows;
            saved.chosen_w = cv_w.chosen;
            saved.train_w = sp_w.train;
            saved.val_w = sp_w.val;
            saved.weights = wmaxtwin_scale_weights(j0, J, weight_mode);
            saved_exists = true;
        end
    end

    records = table(rep_id, Noisy, LevelSure, NasonEvenOddCV, MaxTwinFixedCV, WMaxTwinRampCV, OracleRamp, WMaxTwinSlope, OracleSlope, ...
        'VariableNames', {'rep','Noisy','LevelSure','NasonEvenOddCV','MaxTwinFixedCV','WMaxTwinRampCV','OracleRamp','WMaxTwinSlope','OracleSlope'});

    method_names = {'Noisy'; 'LevelSure'; 'NasonEvenOddCV'; 'MaxTwinFixedCV'; 'WMaxTwinRampCV'; 'OracleRamp'};
    amse = zeros(length(method_names), 1);
    sdv = zeros(length(method_names), 1);
    sev = zeros(length(method_names), 1);
    for ii = 1:length(method_names)
        vals = records.(method_names{ii});
        amse(ii) = mean(vals);
        sdv(ii) = std(vals, 0);
        sev(ii) = sdv(ii) / sqrt(nrep);
    end
    summary = table(method_names, amse, sdv, sev, 'VariableNames', {'Method','AMSE','SD','SE'});

    writetable(records, fullfile(output_dir, 'doppler_replicates_MATLAB.csv'));
    writetable(summary, fullfile(output_dir, 'doppler_amse_summary_MATLAB.csv'));

    save_plots(saved, summary, output_dir, j0, J);

    result.summary = summary;
    result.records = records;
    result.saved = saved;
end

function save_plots(saved, summary, output_dir, j0, J)
    fig = figure('Name', 'One Doppler replication');
    plot(saved.t, saved.y, 'Color', [0.75 0.75 0.75], 'LineWidth', 0.7); hold on;
    plot(saved.t, saved.f, 'k', 'LineWidth', 2.0);
    plot(saved.t, saved.f_sure, 'LineWidth', 1.2, 'LineStyle', '--');
    plot(saved.t, saved.f_nason, 'LineWidth', 1.2, 'LineStyle', ':');
    plot(saved.t, saved.f_maxtwin, 'LineWidth', 1.2, 'LineStyle', '-.');
    plot(saved.t, saved.f_wmax, 'LineWidth', 1.8);
    xlabel('t'); ylabel('signal');
    title('One Doppler replication and fitted reconstructions');
    legend({'noisy','truth','LevelSure','Nason CV','MaxTwin fixed','WMaxTwin ramp'}, 'Location', 'best');
    save_figure(fig, fullfile(output_dir, 'fig03_one_replication.png'));

    fig = figure('Name', 'WMaxTwin split mask');
    train_flag = zeros(length(saved.t), 1);
    train_flag(saved.train_w) = 1;
    scatter(saved.t, train_flag, 12, 'filled');
    xlabel('t'); ylabel('training indicator');
    title('WMaxTwin training-validation mask, first replication');
    ylim([-0.1, 1.1]);
    save_figure(fig, fullfile(output_dir, 'fig04_split_mask.png'));

    fig = figure('Name', 'Selected thresholds');
    js = (j0:(J - 1))';
    plot(js, saved.l_wmax, '-o', 'LineWidth', 1.5); hold on;
    plot(js, saved.l_sure, '--s', 'LineWidth', 1.5);
    xlabel('level j'); ylabel('\lambda_j');
    title('Selected thresholds in the first replication');
    legend({'WMaxTwin ramp','Levelwise SURE'}, 'Location', 'best');
    save_figure(fig, fullfile(output_dir, 'fig05_thresholds.png'));

    fig = figure('Name', 'Ramp-slope CV curve');
    rows = sortrows(saved.rows_w, 'slope');
    plot(rows.slope, rows.score, '-o', 'LineWidth', 1.5); hold on;
    yl = ylim;
    plot([saved.chosen_w.slope(1), saved.chosen_w.slope(1)], yl, '--k', 'LineWidth', 1.0);
    ylim(yl);
    xlabel('ramp slope s'); ylabel('symmetric CV score');
    title('Ramp-slope cross-validation, first replication');
    save_figure(fig, fullfile(output_dir, 'fig06_cv_curve.png'));

    fig = figure('Name', 'AMSE comparison');
    bar(summary.AMSE);
    set(gca, 'XTick', 1:height(summary), 'XTickLabel', summary.Method, 'XTickLabelRotation', 35);
    ylabel('AMSE');
    title('Monte Carlo AMSE comparison');
    save_figure(fig, fullfile(output_dir, 'fig07_amse_bar.png'));

    fig = figure('Name', 'Ramp slope diagnostics');
    plot(1:length(saved.rows_w.slope), saved.rows_w.slope, 'o');
    xlabel('ranked CV candidate'); ylabel('slope s');
    title('Ramp slope candidates, first replication');
    save_figure(fig, fullfile(output_dir, 'fig08_slope_diagnostics.png'));
end

function save_figure(fig, filename)
    set(fig, 'Color', 'w');
    try
        exportgraphics(fig, filename, 'Resolution', 200);
    catch
        print(fig, filename, '-dpng', '-r200');
    end
end
