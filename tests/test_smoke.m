function tests = test_smoke
% Smoke tests for fitcopulatoolbox_allinone (no GUI)
tests = functiontests(localfunctions);
end

function setupOnce(t)
    % --- Robust path bootstrap ---
    restoredefaultpath; rehash toolboxcache;
    repoRoot = fileparts(mfilename('fullpath'));  % folder of this test file
    addpath(genpath(repoRoot));                   % add entire repo (subfolders too)
    clear functions; rehash;

    % Sanity: make sure the main entry exists on path
    assert(exist('fitcopulatoolbox_allinone','file')==2, ...
        'fitcopulatoolbox_allinone.m not found on MATLAB path. ' + ...
        "Ensure the file name matches the function name and it's inside repoRoot or subfolder.");

    % Tiny synthetic example
    rng(42,'twister');
    n = 300;
    x = cumsum(randn(n,1));
    y = 0.7*x + 0.3*randn(n,1);
    t.TestData.T = table(x,y,'VariableNames',{'X1','X2'});
end

function test_end_to_end(t)
    T = t.TestData.T;
    X = T{:,1:2};

    % Get API (nargout>0 -> no GUI)
    api = fitcopulatoolbox_allinone();
    mustHaveFields(api, {'makeUfromMarginals','cluster_copula_module_inline','forecastTauGAS_core'});

    % 1) Marginals -> U
    U = api.makeUfromMarginals(X,'Normal');
    assert(all(isfinite(U(:))) && all(U(:)>0 & U(:)<1), 'U not in (0,1).');

    % 2) Quick copula sanity fit (Gaussian)
    R = copulafit('Gaussian',U);
    assert(issymmetric(R) && all(isfinite(R(:))), 'Bad R from copulafit(Gaussian).');

    % 3) Clustering + per-cluster copula fit (K=2, U-space)
    CC = api.cluster_copula_module_inline(X, ...
        'Marginals','Normal', 'Algorithm','kmeans', ...
        'K',2, 'UseUspace',true, 'Verbose', false);
    assert(isfield(CC,'clusters') && numel(CC.clusters)>=1, 'No clusters returned.');
    assert(all(isfield(CC.clusters(1),{'criteria','fit'})), 'Cluster(1) missing fields.');
    assert(isfinite(double(CC.clusters(1).criteria.BIC)), 'Cluster(1) BIC not finite.');

    % 4) GAS τ forecast
    out = api.forecastTauGAS_core(X, 'Normal', 30, 'Gaussian', 0.10, 0.20);
    m = ~isnan(out.tau_pred);
    assert(any(m) && all(isfinite(out.tau_pred(m))), 'Non-finite tau predictions.');
end

function mustHaveFields(s, names)
    for i=1:numel(names)
        assert(isfield(s, names{i}), 'Missing field: %s', names{i});
    end
end
