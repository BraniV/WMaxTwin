%% WMaxTwin for the Marron-Wand claw density
% This MATLAB script reproduces Example I from the WMaxTwin paper.
%
% The experiment compares three balanced half-sample split geometries:
%
%   1. Random splitting.
%   2. MaxTwin+ splitting based on strengthened nonwavelet local-rank features.
%   3. WMaxTwin+ splitting based on the same MaxTwin+ features augmented with
%      Haar wavelet atoms evaluated at the sample locations.
%
% The density estimator, bandwidth, AMSE loss, data-generating mechanism, and
% observed samples are kept fixed. Only the split geometry changes.
%
% The code is intentionally self-contained. It uses only base MATLAB functions.
% No Wavelet Toolbox and no Statistics and Machine Learning Toolbox are needed.
%
% To use this as a MATLAB notebook, open this file in MATLAB. The %% headings
% define executable sections. In MATLAB Live Editor, choose
%   Save As > Live Script
% to convert it to a native .mlx file.

clear; close all; clc;

%% Reproducibility and output folder
REPS = 150;
SEED = 20260619;
H = 0.11;
OUTPUT_DIR = fullfile(pwd, 'wmaxtwin_claw_matlab_output');

if ~exist(OUTPUT_DIR, 'dir')
    mkdir(OUTPUT_DIR);
end

%% Run the simulation
% The default run uses the same four settings as the paper.
% Because MATLAB and Python have different random-number streams, the exact
% Monte Carlo table will not be bitwise identical to the Python table.
% The expected ordering is the same:
%
%     AMSE(Random) > AMSE(MaxTwin+) > AMSE(WMaxTwin+).

[summary_tbl, reps_tbl, example] = run_simulation(REPS, SEED, H, OUTPUT_DIR);

disp(summary_tbl);

%% Inspect the first few replicate-level results
head(reps_tbl)

%% Output files
% The simulation writes the following files to OUTPUT_DIR:
%
%   wmaxtwin_plus_amse_table_matlab.csv
%   wmaxtwin_plus_replicates_matlab.csv
%   fig_wmaxtwin_plus_split_geometry_matlab.png
%   fig_wmaxtwin_plus_density_estimates_matlab.png
%   fig_wmaxtwin_plus_amse_bars_matlab.png
%   fig_wmaxtwin_plus_amse_differences_matlab.png

fprintf('\nOutput written to:\n  %s\n', OUTPUT_DIR);

%% Main interpretation
% The important point is not that WMaxTwin changes the density estimator. It
% does not. The same half-sample Gaussian KDE and the same bandwidth h = 0.11
% are used for all split geometries. The AMSE changes because the two halves
% are chosen differently.
%
% MaxTwin+ improves over Random by enforcing local rank balance. WMaxTwin+
% improves further by orienting the local order-statistic pairs in a feature
% space that includes Haar scale-location atoms. Those Haar features detect
% narrow local components of the claw density that rank-local geometry alone
% may not fully represent.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Local functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [summary_tbl, reps_tbl, example] = run_simulation(reps, seed, h, output_dir)
    rng(seed, 'twister');

    settings = struct( ...
        'n', {60, 80, 80, 120}, ...
        'J', {5, 5, 6, 5}, ...
        'R', {8, 10, 10, 12});

    grid = linspace(-3.5, 3.5, 701)';
    true_pdf = claw_pdf(grid);

    nrows = reps * numel(settings);
    records = zeros(nrows, 11);
    row = 0;
    example = struct();
    have_example = false;

    for ss = 1:numel(settings)
        n = settings(ss).n;
        J = settings(ss).J;
        R = settings(ss).R;

        for rep = 1:reps
            x = sample_claw(n);

            y_random = random_split(n);
            y_max = maxtwin_plus_split(x, R, 700);
            y_wmax = wmaxtwin_plus_split(x, R, J, 10);

            ar = split_amse(x, y_random, grid, true_pdf, h);
            am = split_amse(x, y_max, grid, true_pdf, h);
            aw = split_amse(x, y_wmax, grid, true_pdf, h);

            row = row + 1;
            records(row, :) = [n, J, R, rep, ar, am, aw, ar - am, am - aw, abs(y_max' * y_wmax) / n, aw - ar];

            if ~have_example && n == 80 && J == 6 && rep == 1
                example.x = x;
                example.y_random = y_random;
                example.y_max = y_max;
                example.y_wmax = y_wmax;
                example.grid = grid;
                example.true_pdf = true_pdf;
                example.h = h;
                have_example = true;
            end
        end
    end

    reps_tbl = array2table(records, 'VariableNames', { ...
        'n', 'J', 'R', 'rep', ...
        'AMSE_Random', 'AMSE_MaxTwin_plus', 'AMSE_WMaxTwin_plus', ...
        'diff_Random_minus_MaxTwin_plus', ...
        'diff_MaxTwin_plus_minus_WMaxTwin_plus', ...
        'abs_signed_overlap', ...
        'diff_WMaxTwin_plus_minus_Random'});

    summary = zeros(numel(settings), 12);
    for ss = 1:numel(settings)
        n = settings(ss).n;
        J = settings(ss).J;
        R = settings(ss).R;
        idx = reps_tbl.n == n & reps_tbl.J == J & reps_tbl.R == R;
        g = reps_tbl(idx, :);

        mr = mean(g.AMSE_Random);
        mm = mean(g.AMSE_MaxTwin_plus);
        mw = mean(g.AMSE_WMaxTwin_plus);
        summary(ss, :) = [ ...
            n, J, R, mr, mm, mw, ...
            100 * (mr - mw) / mr, ...
            100 * (mm - mw) / mm, ...
            100 * (mr - mm) / mr, ...
            median(g.abs_signed_overlap), ...
            mean(g.AMSE_MaxTwin_plus < g.AMSE_Random), ...
            mean(g.AMSE_WMaxTwin_plus < g.AMSE_MaxTwin_plus)];
    end

    summary_tbl = array2table(summary, 'VariableNames', { ...
        'n', 'J', 'R', 'Random', 'MaxTwin_plus', 'WMaxTwin_plus', ...
        'Gain_WMaxTwin_vs_Random_percent', ...
        'Gain_WMaxTwin_vs_MaxTwin_percent', ...
        'Gain_MaxTwin_vs_Random_percent', ...
        'Median_abs_signed_overlap', ...
        'Frac_MaxTwin_better_than_Random', ...
        'Frac_WMaxTwin_better_than_MaxTwin'});

    writetable(summary_tbl, fullfile(output_dir, 'wmaxtwin_plus_amse_table_matlab.csv'));
    writetable(reps_tbl, fullfile(output_dir, 'wmaxtwin_plus_replicates_matlab.csv'));

    make_figures(summary_tbl, reps_tbl, example, output_dir);
end

function y = claw_pdf(x)
    x = x(:);
    y = 0.5 * normal_pdf(x, 0.0, 1.0);
    for ell = 0:4
        y = y + 0.1 * normal_pdf(x, ell / 2.0 - 1.0, 0.1);
    end
end

function x = sample_claw(n)
    probs = [0.5, 0.1, 0.1, 0.1, 0.1, 0.1];
    means = [0.0; -1.0; -0.5; 0.0; 0.5; 1.0];
    sds = [1.0; 0.1; 0.1; 0.1; 0.1; 0.1];
    u = rand(n, 1);
    cp = cumsum(probs);
    comp = zeros(n, 1);
    for i = 1:n
        comp(i) = find(u(i) <= cp, 1, 'first');
    end
    x = means(comp) + sds(comp) .* randn(n, 1);
end

function y = normal_pdf(x, mu, sigma)
    y = exp(-0.5 * ((x - mu) ./ sigma).^2) ./ (sqrt(2 * pi) * sigma);
end

function Fstd = standardize_columns(F)
    F = double(F);
    F = F - mean(F, 1);
    sd = std(F, 0, 1);
    keep = sd > 1.0e-12;
    Fstd = F(:, keep) ./ sd(keep);
end

function u = rank_coordinate(x)
    n = numel(x);
    [~, ord] = sort(x(:));
    ranks = zeros(n, 1);
    ranks(ord) = (1:n)';
    u = ranks / (n + 1.0);
end

function F = maxtwin_plus_features(x, R, M)
    if nargin < 3 || isempty(M)
        M = R;
    end

    x = x(:);
    z = (x - mean(x)) ./ (std(x) + 1.0e-12);
    u = rank_coordinate(x);

    cols = [z, z.^2, z.^3, z.^4];

    for r = 1:R
        lo = (r - 1) / R;
        hi = r / R;
        cols = [cols, double(u > lo & u <= hi)]; %#ok<AGROW>
    end

    tau = 1.2 / (M + 1);
    for m = 1:M
        c = m / (M + 1);
        cols = [cols, exp(-0.5 * ((u - c) / tau).^2)]; %#ok<AGROW>
    end

    F = standardize_columns(cols);
end

function F = haar_wavelet_features(x, J, J0, a, b)
    if nargin < 3 || isempty(J0)
        J0 = max(1, J - 3);
    end
    if nargin < 4 || isempty(a)
        a = -3.0;
    end
    if nargin < 5 || isempty(b)
        b = 3.0;
    end

    x = x(:);
    u = (x - a) / (b - a);
    u = min(max(u, 0.0), 1.0 - 1.0e-12);

    cols = [];
    for j = J0:J
        m = 2^j;
        bins = floor(u * m);
        loc = u * m - bins;
        scale = 2^(j / 2);

        for k = 0:(m - 1)
            v = zeros(size(u));
            idx = bins == k;
            signs = ones(sum(idx), 1);
            signs(loc(idx) >= 0.5) = -1;
            v(idx) = scale * signs;
            cols = [cols, v]; %#ok<AGROW>
        end
    end

    m = 2^J;
    bins = floor(u * m);
    scale = 2^(J / 2);
    for k = 0:(m - 1)
        cols = [cols, scale * double(bins == k)]; %#ok<AGROW>
    end

    F = standardize_columns(cols);
end

function F = wmaxtwin_plus_features(x, R, J)
    Fm = maxtwin_plus_features(x, R, R);
    Fw = haar_wavelet_features(x, J, [], -3.0, 3.0);

    Fm = Fm / sqrt(size(Fm, 2));
    Fw = Fw / sqrt(size(Fw, 2));
    F = [Fm, Fw];
end

function y = random_split(n)
    y = [ones(n / 2, 1); -ones(n / 2, 1)];
    y = y(randperm(n));
end

function y = maxtwin_plus_split(x, R, n_flips)
    x = x(:);
    n = numel(x);
    [~, idx] = sort(x);

    blocks = array_split_indices(idx, R);
    block_id = zeros(n, 1);
    y = zeros(n, 1);
    leftovers = [];

    for bi = 1:R
        bidx = blocks{bi};
        block_id(bidx) = bi;
        bcopy = bidx(randperm(numel(bidx)));
        m = floor(numel(bcopy) / 2);
        y(bcopy(1:m)) = 1;
        y(bcopy((m + 1):(2 * m))) = -1;
        if mod(numel(bcopy), 2) == 1
            leftovers = [leftovers; bcopy(end)]; %#ok<AGROW>
        end
    end

    if ~isempty(leftovers)
        leftovers = leftovers(randperm(numel(leftovers)));
        n_plus = sum(y == 1);
        target_plus = n / 2;
        for ii = 1:numel(leftovers)
            if n_plus < target_plus
                y(leftovers(ii)) = 1;
                n_plus = n_plus + 1;
            else
                y(leftovers(ii)) = -1;
            end
        end
    end

    F = maxtwin_plus_features(x, R, R);
    s = F' * y;
    cur = (s' * s) / size(F, 2);

    for iter = 1:n_flips
        bi = randi(R);
        inds = find(block_id == bi);
        pp = inds(y(inds) == 1);
        mm = inds(y(inds) == -1);
        if isempty(pp) || isempty(mm)
            continue;
        end
        ip = pp(randi(numel(pp)));
        im = mm(randi(numel(mm)));

        ds = -2.0 * F(ip, :)' + 2.0 * F(im, :)';
        ns = s + ds;
        new = (ns' * ns) / size(F, 2);
        if new < cur
            y(ip) = -1;
            y(im) = 1;
            s = ns;
            cur = new;
        end
    end
end

function blocks = array_split_indices(idx, R)
    n = numel(idx);
    q = floor(n / R);
    remn = mod(n, R);
    blocks = cell(R, 1);
    pos = 1;
    for r = 1:R
        sz = q + double(r <= remn);
        blocks{r} = idx(pos:(pos + sz - 1));
        pos = pos + sz;
    end
end

function sgn = orient_adjacent_pairs(pairs, F, restarts, sweeps)
    D = F(pairs(:, 1), :) - F(pairs(:, 2), :);
    m = size(D, 1);
    best_obj = Inf;
    best_sgn = ones(m, 1);

    for rr = 1:restarts
        sgn = 2 * double(rand(m, 1) > 0.5) - 1;
        vec = D' * sgn;
        cur = (vec' * vec) / size(F, 2);

        for sw = 1:sweeps
            improved = false;
            order = randperm(m);
            for ii = 1:m
                q = order(ii);
                new_vec = vec - 2.0 * sgn(q) * D(q, :)';
                new_obj = (new_vec' * new_vec) / size(F, 2);
                if new_obj < cur
                    vec = new_vec;
                    cur = new_obj;
                    sgn(q) = -sgn(q);
                    improved = true;
                end
            end
            if ~improved
                break;
            end
        end

        if cur < best_obj
            best_obj = cur;
            best_sgn = sgn;
        end
    end

    sgn = best_sgn;
end

function y = wmaxtwin_plus_split(x, R, J, restarts)
    x = x(:);
    n = numel(x);
    [~, idx] = sort(x);
    pairs = reshape(idx, 2, [])';

    F = wmaxtwin_plus_features(x, R, J);
    sgn = orient_adjacent_pairs(pairs, F, restarts, 5);

    y = zeros(n, 1);
    for q = 1:size(pairs, 1)
        a = pairs(q, 1);
        b = pairs(q, 2);
        y(a) = sgn(q);
        y(b) = -sgn(q);
    end
end

function f = kde_eval(x, grid, h)
    x = x(:)';
    grid = grid(:);
    Z = (grid - x) / h;
    f = mean(exp(-0.5 * Z.^2) / sqrt(2 * pi), 2) / h;
end

function val = split_amse(x, y, grid, true_pdf, h)
    f_plus = kde_eval(x(y == 1), grid, h);
    f_minus = kde_eval(x(y == -1), grid, h);
    ise_plus = trapz(grid, (f_plus - true_pdf).^2);
    ise_minus = trapz(grid, (f_minus - true_pdf).^2);
    val = 0.5 * (ise_plus + ise_minus);
end

function make_figures(summary_tbl, reps_tbl, example, output_dir)
    x = example.x;
    grid = example.grid;
    true_pdf = example.true_pdf;
    h = example.h;

    %% Figure 1: split geometry
    fig = figure('Color', 'w', 'Position', [100, 100, 850, 400]);
    hold on;
    scatter(x, 0 + 0.12 * example.y_random, 24, 'filled', 'MarkerFaceAlpha', 0.75);
    scatter(x, 1 + 0.12 * example.y_max, 24, 'filled', 'MarkerFaceAlpha', 0.75);
    scatter(x, 2 + 0.12 * example.y_wmax, 24, 'filled', 'MarkerFaceAlpha', 0.75);
    yticks([0, 1, 2]);
    yticklabels({'Random', 'MaxTwin+', 'WMaxTwin+'});
    xlabel('x');
    title('Representative split geometry, n=80, J=6');
    legend({'Random', 'MaxTwin+', 'WMaxTwin+'}, 'Location', 'northeast');
    box on;
    print(fig, fullfile(output_dir, 'fig_wmaxtwin_plus_split_geometry_matlab.png'), '-dpng', '-r200');

    %% Figure 2: half-sample KDE estimates
    fig = figure('Color', 'w', 'Position', [100, 100, 850, 440]);
    hold on;
    plot(grid, true_pdf, 'LineWidth', 2.0);
    plot(grid, kde_eval(x(example.y_max == 1), grid, h), 'LineWidth', 1.1);
    plot(grid, kde_eval(x(example.y_max == -1), grid, h), 'LineWidth', 1.1);
    plot(grid, kde_eval(x(example.y_wmax == 1), grid, h), 'LineWidth', 1.1);
    plot(grid, kde_eval(x(example.y_wmax == -1), grid, h), 'LineWidth', 1.1);
    xlabel('x'); ylabel('density');
    title('Half-sample KDE estimates, n=80, J=6');
    legend({'true density', 'MaxTwin+ A', 'MaxTwin+ B', 'WMaxTwin+ A', 'WMaxTwin+ B'}, ...
        'Location', 'northeast');
    box on;
    print(fig, fullfile(output_dir, 'fig_wmaxtwin_plus_density_estimates_matlab.png'), '-dpng', '-r200');

    %% Figure 3: mean AMSE bar chart
    fig = figure('Color', 'w', 'Position', [100, 100, 880, 440]);
    Y = [summary_tbl.Random, summary_tbl.MaxTwin_plus, summary_tbl.WMaxTwin_plus];
    bar(Y);
    xticklabels(compose('n=%d, J=%d', summary_tbl.n, summary_tbl.J));
    ylabel('mean AMSE');
    title('AMSE comparison across split geometries');
    legend({'Random', 'MaxTwin+', 'WMaxTwin+'}, 'Location', 'northeast');
    box on;
    print(fig, fullfile(output_dir, 'fig_wmaxtwin_plus_amse_bars_matlab.png'), '-dpng', '-r200');

    %% Figure 4: replicatewise AMSE differences
    fig = figure('Color', 'w', 'Position', [100, 100, 880, 440]);
    hold on;
    xs = [];
    labels = {};
    xpos = 1;
    for rr = 1:height(summary_tbl)
        n = summary_tbl.n(rr);
        J = summary_tbl.J(rr);
        R = summary_tbl.R(rr);
        idx = reps_tbl.n == n & reps_tbl.J == J & reps_tbl.R == R;
        dRM = reps_tbl.diff_Random_minus_MaxTwin_plus(idx);
        dMW = reps_tbl.diff_MaxTwin_plus_minus_WMaxTwin_plus(idx);
        draw_simple_box(dRM, xpos);
        labels{end+1} = sprintf('R-M\nn=%d,J=%d', n, J); %#ok<AGROW>
        xs(end+1) = xpos; %#ok<AGROW>
        xpos = xpos + 1;
        draw_simple_box(dMW, xpos);
        labels{end+1} = sprintf('M-W\nn=%d,J=%d', n, J); %#ok<AGROW>
        xs(end+1) = xpos; %#ok<AGROW>
        xpos = xpos + 1;
    end
    yline(0, '--');
    xlim([0.5, xpos - 0.5]);
    xticks(xs);
    xticklabels(labels);
    ylabel('AMSE difference');
    title('Replicatewise AMSE differences: positive values favor the second method');
    box on;
    print(fig, fullfile(output_dir, 'fig_wmaxtwin_plus_amse_differences_matlab.png'), '-dpng', '-r200');
end

function draw_simple_box(data, xpos)
    data = sort(data(:));
    q1 = emp_quantile(data, 0.25);
    q2 = emp_quantile(data, 0.50);
    q3 = emp_quantile(data, 0.75);
    iqr = q3 - q1;
    lo = max(min(data), q1 - 1.5 * iqr);
    hi = min(max(data), q3 + 1.5 * iqr);
    w = 0.30;

    patch([xpos - w, xpos + w, xpos + w, xpos - w], [q1, q1, q3, q3], ...
        [0.85, 0.85, 0.85], 'EdgeColor', 'k');
    line([xpos - w, xpos + w], [q2, q2], 'Color', 'k', 'LineWidth', 1.4);
    line([xpos, xpos], [lo, q1], 'Color', 'k');
    line([xpos, xpos], [q3, hi], 'Color', 'k');
    line([xpos - w/2, xpos + w/2], [lo, lo], 'Color', 'k');
    line([xpos - w/2, xpos + w/2], [hi, hi], 'Color', 'k');
end

function q = emp_quantile(x, p)
    x = sort(x(:));
    n = numel(x);
    if n == 1
        q = x(1);
        return;
    end
    pos = 1 + (n - 1) * p;
    lo = floor(pos);
    hi = ceil(pos);
    if lo == hi
        q = x(lo);
    else
        q = x(lo) + (pos - lo) * (x(hi) - x(lo));
    end
end
