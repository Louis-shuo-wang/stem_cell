function [params_out, info] = calibrate_k1_lambdaP(params_in, idx, targetW, targetRatio)
% ... (前面注释保持不变)

params = params_in;

% ---- user-tunable options (sensible defaults) ----
maxOuterIter = 40;
tolW = 1e-3;
tolR = 1e-3;
calib_Tfinal = max(80, round(params.Tfinal*0.15));
calib_dt = params.dt;
calib_Nt = round(calib_Tfinal / calib_dt);

% override params for calibration short runs
orig_Tfinal = params.Tfinal;
orig_Nt = params.Nt;
params.Tfinal = calib_Tfinal;
params.Nt = calib_Nt;

% k1 bounds (log-space safe)
k1_min = 1e-12;
k1_max = 1e6;

% lambdaP bounds
lambdaP_min = 1e-6;
lambdaP_max = 20;

% initialize
lambdaP = max(lambdaP_min, min(lambdaP_max, params.lambdaP));
k1 = max(k1_min, params.k1);
keep_k2_ratio = true;
k2_ratio = params.k2 / max(params.k1,1e-20);

% bookkeeping
info.iter = 0;
info.history = [];
info.success = false;

% initialize W_now/ratio_now so they always exist
W_now = NaN;
ratio_now = NaN;

% helper: run short sim and return Wsum and ratio
    function [Wsum_out, ratio_out] = run_once_local(p)
        [Psim, Wsim] = main(p, idx);
        Wsum_out = sum(Wsim(:, end)) * p.dx;
        Psum_out = sum(Psim(:, end)) * p.dx;
        ratio_out = Wsum_out / (Psum_out + eps);
    end

% Outer loop: adjust lambdaP
prev_lambdaP = NaN; prev_ratio = NaN;
for outer = 1:maxOuterIter
    info.iter = outer;

    % Inner: find k1 (for current lambdaP) so that W ~= targetW via bisection.
    p_try = params;
    p_try.lambdaP = lambdaP;

    % bracket for k1
    low = k1_min;
    high = max(low*10, k1);
    p_try.k1 = low; if keep_k2_ratio, p_try.k2 = k2_ratio * p_try.k1; end
    [W_low, ~] = run_once_local(p_try);
    p_try.k1 = high; if keep_k2_ratio, p_try.k2 = k2_ratio * p_try.k1; end
    [W_high, ~] = run_once_local(p_try);

    % expand high until W_high < targetW (since W decreases with k1)
    iter_expand = 0;
    while (W_high > targetW) && (high < k1_max) && (iter_expand < 40)
        high = high * 10;
        p_try.k1 = high; if keep_k2_ratio, p_try.k2 = k2_ratio * p_try.k1; end
        [W_high, ~] = run_once_local(p_try);
        iter_expand = iter_expand + 1;
    end

    % flag whether we abandoned the bisection (so we can skip 'for ib' safely)
    bisection_abandoned = false;

    if W_high > targetW
        % failed to bracket; use high bound
        warning('calibrator:unable_to_bracket', 'Failed to bracket W target in k1 up to k1_max. Using high bound.');
        k1_found = high;
        % we'll compute W_now later from params with this k1_found
    else
        % now do bisection on k1 in [low, high]
        k_lo = low; k_hi = high;
        W_lo = W_low; W_hi = W_high;

        % ensure monotone condition: W_lo >= target >= W_hi
        if ~(W_lo >= targetW && W_hi <= targetW)
            if W_lo <= targetW && W_hi >= targetW
                % swap endpoints
                tmp = k_lo; k_lo = k_hi; k_hi = tmp;
                tmpW = W_lo; W_lo = W_hi; W_hi = tmpW;
            else
                % give up on bisection, pick best of low/high
                if abs(W_low - targetW) < abs(W_high - targetW)
                    k1_found = low;
                else
                    k1_found = high;
                end
                info.warning = 'nonmonotonic_W_k1';
                % set into params and evaluate current W and ratio
                params.k1 = k1_found;
                if keep_k2_ratio, params.k2 = k2_ratio * params.k1; end
                [W_now, ratio_now] = run_once_local(params);
                k1 = k1_found;
                bisection_abandoned = true; % skip bisection below
            end
        end

        % perform bisection only if not abandoned
        if ~bisection_abandoned
            max_bisect = 50;
            k1_found = k_lo; % default fallback
            for ib = 1:max_bisect
                k_mid = sqrt(k_lo * k_hi); % geometric mean
                p_try.k1 = k_mid; if keep_k2_ratio, p_try.k2 = k2_ratio * p_try.k1; end
                [W_mid, ~] = run_once_local(p_try);

                if abs(W_mid - targetW) <= tolW
                    k1_found = k_mid;
                    break;
                end

                if W_mid > targetW
                    k_lo = k_mid;
                else
                    k_hi = k_mid;
                end
                k1_found = k_mid;
            end
        end
    end % end bisection/bracket

    % set found k1 into params
    params.k1 = k1_found;
    if keep_k2_ratio, params.k2 = k2_ratio * params.k1; end
    k1 = params.k1;

    % If we didn't already compute W_now (i.e. didn't abandon earlier), compute it now
    if ~bisection_abandoned
        params.lambdaP = lambdaP;
        [W_now, ratio_now] = run_once_local(params);
    end

    info.history(end+1,:) = [outer, params.k1, params.lambdaP, W_now, ratio_now];
    fprintf('Outer %d: k1=%.3e, lambdaP=%.6f, W=%.6f, ratio=%.6f\n', outer, params.k1, params.lambdaP, W_now, ratio_now);

    % check convergence
    if abs(W_now - targetW) <= tolW && abs(ratio_now - targetRatio) <= tolR
        info.success = true;
        break;
    end

    % Update lambdaP using multiplicative control (robust)
    if ~isnan(ratio_now) && ratio_now > 0
        factor = (targetRatio / (ratio_now + eps)) ^ 0.6;
    else
        factor = 1.2;
    end

    lambdaP_new = lambdaP * factor;
    lambdaP_new = max(lambdaP_min, min(lambdaP_max, lambdaP_new));
    lambdaP_new = max(lambdaP * 0.5, min(lambdaP * 1.5, lambdaP_new));

    prev_lambdaP = lambdaP;
    prev_ratio = ratio_now;
    lambdaP = lambdaP_new;
    params.lambdaP = lambdaP;
end % outer loop

% restore original long-run Tfinal/Nt
params.Tfinal = orig_Tfinal;
params.Nt = orig_Nt;

% ensure final W_now/ratio_now exist (safety)
if isnan(W_now)
    [W_now, ratio_now] = run_once_local(params);
end

params_out = params;
info.finalW = W_now;
info.finalRatio = ratio_now;
info.k1 = params.k1;
info.lambdaP = params.lambdaP;
info.targetW = targetW;
info.targetRatio = targetRatio;
end
