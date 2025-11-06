function tests = test_smoke
tests = functiontests(localfunctions);
end

function setupOnce(t)
% Fix seeds for reproducibility
rng(42,'twister');
t.TestData.example = load_example();        % helper below
t.TestData.cols = t.TestData.example.Properties.VariableNames(1:2);
end

function test_end_to_end(~)
% Load data
T = load_example();                         %#ok<NASGU>

% Marginals -> U
X = T{:,1:2};
U = makeUfromMarginals(X,'Normal');         % uses toolbox function
assert(all(isfinite(U(:))) && all(U(:)>0 & U(:)<1));

% Fit a copula (Gaussian or t)
R = copulafit('Gaussian',U);
assert(issymmetric(R) && all(isfinite(R(:))));

% Quick clustering + per-cluster fit
CC = cluster_copula_module_inline(X, 'Marginals','Normal', ...
    'Algorithm','kmeans','K',2,'UseUspace',true, 'Verbose', false);
assert(isfield(CC,'clusters') && numel(CC.clusters)>=1);
assert(all(isfinite(double(CC.clusters(1).criteria.BIC))));

% GAS τ forecast on a short rolling window
out = local_forecast_tau_gas(X, 'Normal', 30); % helper below
m = ~isnan(out.tau_pred);
assert(any(m) && all(isfinite(out.tau_pred(m))));
end

%% ----------------- local helpers for the test -----------------
function T = load_example()
% Use a tiny, self-contained example (or read from examples/demo.csv)
n = 300; x = cumsum(randn(n,1)); y = 0.7*x + 0.3*randn(n,1);
T = table(x,y,'VariableNames',{'X1','X2'});
end

function out = local_forecast_tau_gas(X, marg, W)
% Thin wrapper to call the toolbox’s GAS core with a tiny window
out = forecastTauGAS_core(X, marg, W, 'Gaussian', 0.10, 0.2);
end
