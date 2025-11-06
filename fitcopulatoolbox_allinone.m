
function api = fitcopulatoolbox_allinone(varargin)
% FITCOPULATOOLBOX_ALLINONE — Single-file Copula Toolbox with Clustering & Mixtures
% -----------------------------------------------------------------------------
% GUI για:
%  - Φόρτωση δεδομένων
%  - Επιλογή marginals και fit Copula (Gaussian / t / Clayton / Frank / Gumbel)
%  - Διαγνωστικά (PIT, heatmaps, κ.λπ.)
%  - Cluster analysis (k-means / GMM / Hierarchical / DBSCAN) πάνω σε X ή U-space
%  - Ανά cluster fit copula (model selection με BIC)
%  - Προαιρετικά Mixture-of-copulas (2D) με EM
%  - Εξαγωγές / Save-Load sessions
%  - Συμβατό με MATLAB R2023b+
%
% Σημείωση: Το παρόν αρχείο συνδυάζει ό,τι χρειάζεται σε μία συνάρτηση αρχείου.
% Όλες οι βοηθητικές ρουτίνες ορίζονται ως nested subfunctions.
%
% Συγγραφείς: Statiou D. Anastasios & Hatzopoulos Petros, 2025
% Άδεια χρήσης: Επιτρεπτή ενσωμάτωση/τροποποίηση με αναφορά.
%
% -----------------------------------------------------------------------------

% FITCOPULATOOLBOX_ALLINONE — Single-file Copula Toolbox with Clustering & Mixtures

% ---- API mode ----
if nargout >= 1 || (nargin>=1 && ischar(varargin{1}) && strcmpi(varargin{1},'api'))
    api = struct( ...
        'makeUfromMarginals',           @makeUfromMarginals, ...
        'cluster_copula_module_inline', @cluster_copula_module_inline, ...
        'forecastTauGAS_core',          @forecastTauGAS_core ...
    );
    return;
end
%% Main GUI
f = uifigure('Name','Copula Toolbox — All-in-One', 'Position',[40 40 1360 780]);

% Left panel controls ------------------------------------------------------
uilabel(f,'Position',[20 730 120 22],'Text','Copula:');
ddlCopula = uidropdown(f,'Position',[120 730 230 24], ...
    'Items',{'Gaussian','t','Clayton (2D)','Frank (2D)','Gumbel (2D)'}, 'Value','Gaussian');

uilabel(f,'Position',[20 700 120 22],'Text','Marginals:');
ddlMarg = uidropdown(f,'Position',[120 700 230 24], ...
    'Items',{'Normal','Lognormal','Exponential','Gamma','tLocationScale','Ranks (empirical)'}, ...
    'Value','Normal');

btnLoad = uibutton(f,'push','Text','Load Data','Position',[370 730 120 26], ...
    'ButtonPushedFcn', @(src,~) loadData(src));
uilabel(f,'Position',[20 670 120 22],'Text','Columns:');
lstCols = uilistbox(f,'Position',[120 530 230 130],'Multiselect','on','Items',{},'Value',string.empty(0,1));

% Numeric indices (optional)
uilabel(f,'Position',[20 500 120 22],'Text','Col idx:');
edtIdx  = uieditfield(f,'text','Position',[120 500 230 24],'Placeholder','e.g. 1,2,3','Value','');

% Extra controls
uilabel(f,'Position',[20 470 120 22],'Text','alpha:');
edtAlpha = uieditfield(f,'numeric','Position',[120 470 80 24],'Value',0.10,'Limits',[0.001 0.5]);
uilabel(f,'Position',[210 470 60 22],'Text','K-fold:');
edtKfold = uieditfield(f,'numeric','Position',[270 470 80 24],'Value',5,'Limits',[2 20],'RoundFractionalValues',true);
uilabel(f,'Position',[20 440 120 22],'Text','Bootstrap B:');
edtBootB = uieditfield(f,'numeric','Position',[120 440 80 24],'Value',200,'Limits',[50 5000],'RoundFractionalValues',true);

% Action buttons row 1
uibutton(f,'push','Text','Fit Copula','Position',[510 730 120 26],...
    'ButtonPushedFcn', @(src,~) fitMultivariate(src));
uibutton(f,'push','Text','Scatter Matrix','Position',[640 730 130 26],...
    'ButtonPushedFcn', @(src,~) showScatterMatrix(src));
uibutton(f,'push','Text','Compare R-matrices','Position',[780 730 150 26],...
    'ButtonPushedFcn', @(src,~) compareCorrelations(src));
uibutton(f,'push','Text','Marginal Histograms','Position',[940 730 150 26],...
    'ButtonPushedFcn', @(src,~) showMarginals(src));
uibutton(f,'push','Text','Auto-Select Copula','Position',[1100 730 150 26],...
    'ButtonPushedFcn', @(src,~) autoSelectCopula(src));

% Action buttons row 2
uibutton(f,'push','Text','Joint PDF (2D)','Position',[510 700 120 26],...
    'ButtonPushedFcn', @(src,~) plotJointPDF(src));
uibutton(f,'push','Text','Tail Dependence','Position',[640 700 130 26],...
    'ButtonPushedFcn', @(src,~) tailDependenceDialog(src));
uibutton(f,'push','Text','GOF (Rosenblatt, 2D)','Position',[780 700 150 26],...
    'ButtonPushedFcn', @(src,~) gofDialog(src));
uibutton(f,'push','Text','Rolling tau','Position',[940 700 150 26],...
    'ButtonPushedFcn', @(src,~) rollingTauDialog(src));
uibutton(f,'push','Text','Forecast tau(t+1)','Position',[1100 700 150 26],...
    'ButtonPushedFcn', @(src,~) forecastDialog(src));

% Action buttons row 3
uibutton(f,'push','Text','PIT Diagnostics','Position',[510 670 120 26],...
    'ButtonPushedFcn', @(src,~) pitDiagnostics(src));
uibutton(f,'push','Text','Rank Heatmaps','Position',[640 670 130 26],...
    'ButtonPushedFcn', @(src,~) rankHeatmaps(src));
uibutton(f,'push','Text','Empirical Copula (2D)','Position',[780 670 150 26],...
    'ButtonPushedFcn', @(src,~) empiricalCopulaContours(src));
uibutton(f,'push','Text','Copula PDF on U (2D)','Position',[940 670 150 26],...
    'ButtonPushedFcn', @(src,~) copulaPDFonU(src));
uibutton(f,'push','Text','Copula CDF on U (2D)','Position',[1100 670 150 26],...
    'ButtonPushedFcn', @(src,~) copulaCDFonU(src));

% Action buttons row 4
uibutton(f,'push','Text','Kendall Plot (2D)','Position',[510 640 120 26],...
    'ButtonPushedFcn', @(src,~) kendallPlot2D(src));
uibutton(f,'push','Text','K-fold CV LogLik','Position',[640 640 130 26],...
    'ButtonPushedFcn', @(src,~) kfoldCV(src));
uibutton(f,'push','Text','Param Bootstrap','Position',[780 640 150 26],...
    'ButtonPushedFcn', @(src,~) bootstrapDialog(src));
uibutton(f,'push','Text','3D Scatter (2D/3D)','Position',[940 640 150 26],...
    'ButtonPushedFcn', @(src,~) Uscatter3D(src));
uibutton(f,'push','Text','Conditional Slice (2D)','Position',[1100 640 150 26],...
    'ButtonPushedFcn', @(src,~) conditionalSliceDialog(src));

% Sim/Export row
uilabel(f,'Position',[510 610 40 22],'Text','nSim:');
edtNSim = uieditfield(f,'numeric','Position',[550 610 80 24],'Value',1000,'Limits',[10 Inf],'RoundFractionalValues',true);
uibutton(f,'push','Text','Simulate from Copula','Position',[640 610 170 26],...
    'ButtonPushedFcn', @(src,~) simulateFromCopula(src));
uibutton(f,'push','Text','Export U/R/fit','Position',[820 610 130 26],...
    'ButtonPushedFcn', @(src,~) exportResults(src));
uibutton(f,'push','Text','Save Session','Position',[960 610 120 26],...
    'ButtonPushedFcn', @(src,~) saveSession(src));
uibutton(f,'push','Text','Load Session','Position',[1090 610 120 26],...
    'ButtonPushedFcn', @(src,~) loadSession(src));

% Export graphics / report
uibutton(f,'push','Text','Export Axes PNG','Position',[30 415 160 26],...
    'ButtonPushedFcn', @(src,~) exportAxesPNG(src));
uibutton(f,'push','Text','Quick Report','Position',[190 415 160 26],...
    'ButtonPushedFcn', @(src,~) quickReport(src));

% ---- Clustering controls -------------------------------------------------
uilabel(f,'Position',[20 390 120 22],'Text','Clustering:');
ddlClAlg = uidropdown(f,'Position',[120 390 230 24], ...
    'Items',{'kmeans','gmm','hier','dbscan'}, 'Value','kmeans');

uilabel(f,'Position',[20 360 120 22],'Text','K (blank=auto):');
edtClK = uieditfield(f,'text','Position',[120 360 230 24],'Value','');

uilabel(f,'Position',[20 330 120 22],'Text','Kmax:');
edtClKmax = uieditfield(f,'numeric','Position',[120 330 80 24],'Value',8,'Limits',[2 50],'RoundFractionalValues',true);

uilabel(f,'Position',[210 330 60 22],'Text','Use U:');
chkClUseU = uicheckbox(f,'Position',[270 330 80 24],'Value',false);

uilabel(f,'Position',[20 300 120 22],'Text','DBSCAN ε / MinPts:');
edtClEps = uieditfield(f,'numeric'); edtClEps.Position = [120 300 80 24]; edtClEps.Limits = [0 Inf]; edtClEps.Value = 0;
edtClMin = uieditfield(f,'numeric'); edtClMin.Position = [210 300 80 24]; edtClMin.Limits = [0 Inf]; edtClMin.Value = 0;

uilabel(f,'Position',[20 270 120 22],'Text','Mixture 2D:');
chkMix2D = uicheckbox(f,'Position',[120 270 80 24],'Value',false);
uilabel(f,'Position',[210 270 60 22],'Text','MaxK:');
edtMixK = uieditfield(f,'numeric','Position',[270 270 80 24],'Value',3,'Limits',[2 8],'RoundFractionalValues',true);

% κουμπιά clustering
uilabel(f,'Position',[20 260 240 22],'Text','CLUSTERING — actions','FontWeight','bold');
uibutton(f,'push','Text','Cluster & Fit Copulas','Position',[20 240 160 26], ...
    'ButtonPushedFcn', @(src,~) runClusterAndFit(src));

% Extra Cluster Analysis buttons (added)
uibutton(f,'push','Text','Auto K (grid)','Position',[20 180 160 26],...
    'ButtonPushedFcn', @(src,~) autoKGridDialog(src));
uibutton(f,'push','Text','DBSCAN ε elbow','Position',[190 90 170 26],...
    'ButtonPushedFcn', @(src,~) dbscanElbowDialog(src));
uibutton(f,'push','Text','Outlier Prune (U)','Position',[20 150 160 26],...
    'ButtonPushedFcn', @(src,~) outlierPruneUDialog(src));
uibutton(f,'push','Text','Per-cluster GOF','Position',[190 150 170 26],...
    'ButtonPushedFcn', @(src,~) perClusterGOFDialog(src));
uibutton(f,'push','Text','Feature Importance','Position',[20 120 160 26],...
    'ButtonPushedFcn', @(src,~) featureImportanceDialog(src));
uibutton(f,'push','Text','Export Cluster Report','Position',[190 120 170 26],...
    'ButtonPushedFcn', @(src,~) exportClusterReport(src));
uibutton(f,'push','Text','Spectral (Copula CvM)','Position',[20 90 160 26],...
    'ButtonPushedFcn', @(src,~) spectralCopulaCvMDialog(src));

uibutton(f,'push','Text','Show Cluster Plots','Position',[190 240 160 26], ...
    'ButtonPushedFcn', @(src,~) showClusterPlots(src));
uibutton(f,'push','Text','Export Cluster Results','Position',[20 210 160 26], ...
    'ButtonPushedFcn', @(src,~) exportClusterResults(src));
% extra analysis buttons
uibutton(f,'push','Text','HSIC Independence','Position',[190 210 170 26], ...
    'ButtonPushedFcn', @(src,~) hsicTestDialog(src));
uibutton(f,'push','Text','CoVaR/CoES (2D)','Position',[20 60 160 26], ...
    'ButtonPushedFcn', @(src,~) covarDialog(src));
uibutton(f,'push','Text','Cluster Stability','Position',[190 180 160 26], ...
    'ButtonPushedFcn', @(src,~) clusterStabilityDialog(src));
uibutton(f,'push','Text','Elbow/Silhouette','Position',[190 60 160 26], ...
    'ButtonPushedFcn', @(src,~) elbowSilhouette(src));


% Right: axes + output
ax = uiaxes(f,'Position',[380 140 940 460]); title(ax, 'Correlation Heatmap');
txtOut = uitextarea(f, 'Position',[380 40 940 100], 'Editable','off','FontName','Consolas');
uilabel(f,'Position',[980 740 320 50],'Text','Statiou D. A. & Hatzopoulos P.','FontAngle','italic','FontWeight','bold','HorizontalAlignment','right');

% State
userData = struct('Table',[], 'U',[], 'R',[], 'nu',[], 'theta',[], 'family','', 'marg','Normal', ...
                  'Cluster',struct('results',[], 'lastAlg','', 'lastK',[], 'lastUseU',false));

try, close(splash); end

%% ===================== Helpers (shared) =====================
    function cols = getSelectedColumns()
        cols = [];
        if ~isempty(strtrim(edtIdx.Value)) && ~isempty(userData.Table)
            txt  = regexprep(edtIdx.Value,'\[s, ','');
            nums = regexp(txt,'\d+','match');
            if ~isempty(nums)
                idx = unique(str2double(nums));
                idx = idx(~isnan(idx) & idx>=1 & idx<=width(userData.Table));
                if ~isempty(idx)
                    cols = userData.Table.Properties.VariableNames(idx);
                end
            end
        end
        if isempty(cols), cols = lstCols.Value; end
        if isstring(cols) || iscategorical(cols), cols = cellstr(cols); end
        if ischar(cols), cols = {cols}; end
        cols = cols(:);
    end

    function U = makeUfromMarginals(X, margName)
        n = size(X,1); d = size(X,2);
        U = zeros(n,d);
        if strcmpi(margName,'Ranks (empirical)')
            for j=1:d
                [~,r] = sort(X(:,j)); invr = zeros(n,1); invr(r) = 1:n; U(:,j) = invr/(n+1);
            end
        else
            for j=1:d
                try
                    pd = fitdist(X(:,j), lower(margName));
                    U(:,j) = cdf(pd, X(:,j));
                catch
                    [~,r] = sort(X(:,j)); invr = zeros(n,1); invr(r) = 1:n; U(:,j) = invr/(n+1);
                end
            end
        end
        U = min(max(U, 1e-12), 1-1e-12);
    end

    function assertFitted(src)
        if isempty(userData.U)
            showAlert(src,'Fit a copula first.','Info');
            error('Not fitted.');
        end
    end

%% ===================== Load Data =====================
    function loadData(src)
        try
            [file, path] = uigetfile({'*.csv;*.xlsx;*.xls;*.mat','Data Files (*.csv, *.xlsx, *.xls, *.mat)'}, 'Select Data File');
            if isequal(file,0), return; end
            fullpath = fullfile(path,file); [~,~,ext] = fileparts(fullpath);
            switch lower(ext)
                case '.csv'
                    T = readtable(fullpath,'PreserveVariableNames',true);
                case {'.xlsx','.xls'}
                    T = readtable(fullpath,'PreserveVariableNames',true);
                case '.mat'
                    s = load(fullpath); fn = fieldnames(s); T = s.(fn{1});
                otherwise
                    showAlert(src,'Unsupported file format','Error'); return;
            end
            if ~istable(T), T = array2table(T); end
            names = T.Properties.VariableNames(:);
            if isstring(names) || iscategorical(names), names = cellstr(names); end
            names = cellfun(@char, names, 'UniformOutput', false);
            userData.Table = T; lstCols.Value = string.empty(0,1); lstCols.Items = names;
            if numel(names) >= 2, lstCols.Value = names(1:2); else, lstCols.Value = names(1); end
            edtIdx.Value = '';
            txtOut.Value = ""; userData.U=[]; userData.R=[]; userData.nu=[]; userData.theta=[]; userData.family='';
            cla(ax); title(ax,'Correlation Heatmap');
        catch ME
            showAlert(src, ME.message, 'Load Error');
        end
    end

%% ===================== Fit (main) =====================
    function fitMultivariate(src)
        cols = getSelectedColumns();
        if isempty(cols) || numel(cols) < 2
            showAlert(src,'Select at least 2 columns','Error'); return;
        end
        X = userData.Table{:,cols};
        marg = ddlMarg.Value; userData.marg = marg;
        U = makeUfromMarginals(X, marg);
        copulaType = ddlCopula.Value;
        try
            switch lower(copulaType)
                case 'gaussian'
                    R = copulafit('Gaussian', U);
                    R = nearestPD(R); nu = []; theta = []; famKey = 'Gaussian';
                case 't'
                    [R, nu] = copulafit('t', U);
                    R = nearestPD(R); theta = []; famKey = 't';
                case 'clayton (2d)'
                    if size(U,2) ~= 2, showAlert(src,'Clayton supports only 2D.','Copula Restriction'); return; end
                    theta = copulafit('Clayton', U); R=[]; nu=[]; famKey='Clayton';
                case 'frank (2d)'
                    if size(U,2) ~= 2, showAlert(src,'Frank supports only 2D.','Copula Restriction'); return; end
                    theta = copulafit('Frank', U); R=[]; nu=[]; famKey='Frank';
                case 'gumbel (2d)'
                    if size(U,2) ~= 2, showAlert(src,'Gumbel supports only 2D.','Copula Restriction'); return; end
                    theta = copulafit('Gumbel', U); R=[]; nu=[]; famKey='Gumbel';
                otherwise
                    showAlert(src,'Unsupported copula','Error'); return;
            end
        catch ME
            showAlert(src, ME.message, 'Fit Error'); return;
        end

        userData.U = U; userData.R = R; userData.theta = theta; userData.nu = nu; userData.family = famKey;
        txt = sprintf('Fitted %s copula on %d variables\nMarginals: %s\n', famKey, size(U,2), marg);
        if ~isempty(nu), txt = [txt, sprintf('Degrees of freedom (nu): %.4f\n', nu)]; end
        if ~isempty(R)
            txt = [txt, 'R =', newline, evalc('disp(R)')];
            imagesc(ax, R); colorbar(ax);
            xticks(ax, 1:numel(cols)); yticks(ax, 1:numel(cols)); xticklabels(ax, cols); yticklabels(ax, cols);
            title(ax,'Copula Correlation (R)');
        else
            cla(ax); title(ax, sprintf('%s fitted (2D)', famKey));
        end
        txtOut.Value = txt;
    end

%% ===================== Plots / Analysis =====================
    function showScatterMatrix(src)
        if isempty(userData.Table)
            showAlert(src,'Load data first','Error'); return;
        end
        cols = getSelectedColumns();
        if numel(cols) < 2, showAlert(src,'Select ≥2 columns','Error'); return; end
        X = userData.Table{:,cols};
        figure('Name','Scatter Matrix'); plotmatrix(X); sgtitle('Scatter Matrix of Selected Variables');
    end

    function compareCorrelations(src)
        if isempty(userData.Table)
            showAlert(src,'Load data first','Error'); return;
        end
        if any(strcmpi(userData.family,{'Clayton','Frank','Gumbel'}))
            showAlert(src,'Compare R available only for Gaussian/t.','Info'); return;
        end
        cols = getSelectedColumns();
        if numel(cols)<2, showAlert(src,'Select ≥2 columns','Error'); return; end
        X = userData.Table{:,cols};
        U = makeUfromMarginals(X, ddlMarg.Value);
        try
            switch userData.family
                case 't'
                    [Rsub,~] = copulafit('t', U);
                otherwise
                    Rsub = copulafit('Gaussian', U);
            end
        catch
            Rsub = corr(norminv(U), 'rows','pairwise');
        end
        pearsonR = corr(X, 'type','Pearson','rows','pairwise');
        figure('Name','Correlation Comparison');
        subplot(1,2,1); imagesc(pearsonR); colorbar; title('Pearson Correlation');
        xticks(1:size(X,2)); yticks(1:size(X,2)); xticklabels(cols); yticklabels(cols);
        subplot(1,2,2); imagesc(Rsub); colorbar; title('Copula (R) [subset]');
        xticks(1:size(X,2)); yticks(1:size(X,2)); xticklabels(cols); yticklabels(cols);
    end

    function showMarginals(src)
        if isempty(userData.Table), showAlert(src,'Load data first','Error'); return; end
        cols = getSelectedColumns(); if isempty(cols), showAlert(src,'Select columns','Error'); return; end
        X = userData.Table{:,cols}; marg = ddlMarg.Value;
        figure('Name','Marginals'); m = size(X,2); r = ceil(sqrt(m)); c = ceil(m/r);
        for i = 1:m
            subplot(r,c,i); histogram(X(:,i),'Normalization','pdf'); hold on;
            if ~strcmpi(marg,'Ranks (empirical)')
                try
                    pd = fitdist(X(:,i), lower(marg)); x_range = linspace(min(X(:,i)), max(X(:,i)), 200);
                    plot(x_range, pdf(pd, x_range), 'LineWidth',1.4);
                catch
                end
            end
            title(cols{i},'Interpreter','none');
        end
        sgtitle(sprintf('Marginal Histograms (marginals: %s)', marg));
    end

    function rankHeatmaps(src)
        if isempty(userData.Table), showAlert(src,'Load data first','Error'); return; end
        cols = getSelectedColumns(); if numel(cols)<2, showAlert(src,'Select ≥2','Error'); return; end
        X = userData.Table{:,cols};
        K = corr(X,'type','Kendall','rows','pairwise');
        S = corr(X,'type','Spearman','rows','pairwise');
        figure('Name','Rank Correlations');
        subplot(1,2,1); imagesc(K); colorbar; title('Kendall tau');
        xticks(1:numel(cols)); yticks(1:numel(cols)); xticklabels(cols); yticklabels(cols);
        subplot(1,2,2); imagesc(S); colorbar; title('Spearman \rho');
        xticks(1:numel(cols)); yticks(1:numel(cols)); xticklabels(cols); yticklabels(cols);
    end

    function pitDiagnostics(src)
        if isempty(userData.Table), showAlert(src,'Load data first','Error'); return; end
        cols = getSelectedColumns(); if numel(cols)<1, showAlert(src,'Select columns','Error'); return; end
        X = userData.Table{:,cols}; U = makeUfromMarginals(X, ddlMarg.Value);
        m = size(U,2); r = ceil(sqrt(m)); c = ceil(m/r);
        figure('Name','PIT Diagnostics (U~U[0,1])');
        for j=1:m
            subplot(r,c,j);
            histogram(U(:,j), 'Normalization','pdf'); hold on;
            x = linspace(0,1,200); plot(x, ones(size(x)), 'LineWidth',1.2);
            title(sprintf('%s (PIT)', cols{j}), 'Interpreter','none'); xlim([0 1]);
        end
        sgtitle('PIT histograms: good marginals => approx. flat at 1');
    end

    function empiricalCopulaContours(src)
        if isempty(userData.Table), showAlert(src,'Load data first','Error'); return; end
        cols = getSelectedColumns(); if numel(cols)~=2, showAlert(src,'Select exactly 2 columns','Empirical Copula'); return; end
        X = userData.Table{:,cols}; U = makeUfromMarginals(X, ddlMarg.Value); U1=U(:,1); U2=U(:,2);
        g = linspace(0,1,60); [G1,G2] = meshgrid(g,g);
        Cemp = zeros(size(G1));
        for i=1:numel(G1)
            Cemp(i) = mean(U1<=G1(i) & U2<=G2(i));
        end
        figure('Name','Empirical Copula (CDF)');
        contourf(g,g,Cemp, 12, 'LineStyle','none'); colorbar;
        xlabel('u'); ylabel('v'); title(sprintf('Empirical Copula C_n(u,v): %s-%s', cols{1}, cols{2}));
    end

    function copulaPDFonU(src)
        assertFitted(src);
        cols = getSelectedColumns();
        if numel(cols)~=2, showAlert(src,'Select exactly 2 columns','PDF on U'); return; end
        Xpair = userData.Table{:,cols};
        Upair = makeUfromMarginals(Xpair, ddlMarg.Value);
        fam = userData.family;
        switch fam
            case 't'
                [R2,nu2] = copulafit('t', Upair); param = {R2, max(2,round(nu2))};
            case 'Gaussian'
                R2 = copulafit('Gaussian', Upair); param = {R2};
            otherwise
                th2 = copulafit(fam, Upair); param = {th2};
        end
        g = linspace(0,1,80); [G1,G2] = meshgrid(g,g);
        switch fam
            case 't',        Cpdf = copulapdf('t', [G1(:) G2(:)], param{1}, param{2});
            case 'Gaussian', Cpdf = copulapdf('Gaussian', [G1(:) G2(:)], param{1});
            otherwise,       Cpdf = copulapdf(fam, [G1(:) G2(:)], param{1});
        end
        Cpdf = reshape(Cpdf,size(G1));
        figure('Name','Copula PDF on U'); surf(g,g,Cpdf,'EdgeColor','none'); view(135,30); colorbar;
        xlabel('u'); ylabel('v'); zlabel('c(u,v)');
        title(sprintf('%s copula PDF on U — pair: %s-%s', fam, cols{1}, cols{2}));
    end

    function copulaCDFonU(src)
        assertFitted(src);
        cols = getSelectedColumns(); if numel(cols)~=2, showAlert(src,'Select exactly 2 columns','CDF on U'); return; end
        Xpair = userData.Table{:,cols}; Upair = makeUfromMarginals(Xpair, ddlMarg.Value);
        fam = userData.family;
        switch fam
            case 't'
                [R2,nu2] = copulafit('t', Upair); param = {R2, max(2,round(nu2))};
            case 'Gaussian'
                R2 = copulafit('Gaussian', Upair); param = {R2};
            otherwise
                th2 = copulafit(fam, Upair); param = {th2};
        end
        g = linspace(0,1,80); [G1,G2] = meshgrid(g,g);
        switch fam
            case 't',        Ccdf = copulacdf('t', [G1(:) G2(:)], param{1}, param{2});
            case 'Gaussian', Ccdf = copulacdf('Gaussian', [G1(:) G2(:)], param{1});
            otherwise,       Ccdf = copulacdf(fam, [G1(:) G2(:)], param{1});
        end
        Ccdf = reshape(Ccdf,size(G1));
        figure('Name','Copula CDF on U'); contourf(g,g,Ccdf, 12, 'LineStyle','none'); colorbar;
        xlabel('u'); ylabel('v'); title(sprintf('%s copula CDF on U — pair: %s-%s', fam, cols{1}, cols{2}));
    end

    function kendallPlot2D(src)
        assertFitted(src);
        cols = getSelectedColumns(); if numel(cols)~=2, showAlert(src,'Select exactly 2 columns','Kendall Plot'); return; end
        X = userData.Table{:,cols}; U = makeUfromMarginals(X, ddlMarg.Value);
        fam = userData.family;
        switch fam
            case 't'
                [R2,nu2] = copulafit('t', U); W = copulacdf('t',U, R2, max(2,round(nu2)));
            case 'Gaussian'
                R2 = copulafit('Gaussian', U); W = copulacdf('Gaussian',U, R2);
            otherwise
                th2 = copulafit(fam, U); W = copulacdf(fam, U, th2);
        end
        [femp,xx] = ecdf(W);
        t = linspace(0,1,200); K0 = t - t.*log(max(t,eps));
        figure('Name','Kendall Plot');
        plot(xx,femp,'LineWidth',1.6); hold on; plot(t,K0,'--','LineWidth',1.2); grid on;
        xlabel('t'); ylabel('K(t)'); legend('Empirical K','Independence K_0(t)=t-t\log t','Location','best');
        title(sprintf('Kendall plot: %s-%s (%s)', cols{1}, cols{2}, fam));
    end

%% ===================== Tail dependence =====================
    function tailDependenceDialog(src)
        cols = getSelectedColumns(); if numel(cols) < 2, showAlert(src,'Select ≥2 columns','Tail'); return; end
        if isempty(userData.family), showAlert(src,'Fit a copula first.','Tail'); return; end
        fam = userData.family;
        if any(strcmpi(fam,{'Clayton','Frank','Gumbel'})) && numel(cols)~=2
            showAlert(src,'Tail dep. for Archimedeans is 2D-only here. Select 2 columns.','Tail'); return;
        end
        if any(strcmpi(fam,{'Gaussian','t'}))
            X = userData.Table{:,cols};
            U = makeUfromMarginals(X, ddlMarg.Value);
            try
                if strcmpi(fam,'t')
                    [Rsub,nu2] = copulafit('t', U);
                else
                    Rsub = copulafit('Gaussian', U);
                end
            catch
                Rsub = corr(norminv(U),'rows','pairwise');
            end
            names = cols; d = numel(cols); pairs = nchoosek(1:d,2);
            out = cell(size(pairs,1), 4);
            for i=1:size(pairs,1)
                a = pairs(i,1); b = pairs(i,2); rho = Rsub(a,b);
                if strcmpi(fam,'Gaussian')
                    lamL=0; lamU=0;
                else
                    nuv = max(2,round(nu2));
                    tval = -sqrt((nuv+1)*(1-rho)/(1+rho));
                    lam = 2*tcdf(tval, nuv+1);
                    lamL = lam; lamU = lam;
                end
                out(i,:) = {sprintf('%s-%s',names{a},names{b}), rho, lamL, lamU};
            end
        else
            theta = userData.theta; switch fam
                case 'Clayton', lamL = 2^(-1/theta); lamU = 0;
                case 'Gumbel',  lamL = 0; lamU = 2 - 2^(1/theta);
                otherwise,      lamL = 0; lamU = 0; % Frank
            end
            out = {sprintf('%s-%s',cols{1},cols{2}), NaN, lamL, lamU};
        end
        fig = uifigure('Name','Tail Dependence','Position',[360 360 560 260]);
        uitable(fig,'Data',out,'ColumnName',{'Pair','rho','lambda_L','lambda_U'},'Position',[20 20 520 210]);
    end

%% ===================== GOF (2D Rosenblatt) =====================
    function gofDialog(src)
        cols = getSelectedColumns();
        if numel(cols)~=2, showAlert(src,'Select exactly 2 columns.','GOF'); return; end
        if isempty(userData.family), showAlert(src,'Fit a copula first.','GOF'); return; end
        runGOF(src, cols);
    end

    function runGOF(src, cols)
        X = userData.Table{:,cols}; marg = ddlMarg.Value; U = makeUfromMarginals(X, marg);
        fam = userData.family; epsFD = 1e-4;
        switch fam
            case 'Gaussian'
                R = copulafit('Gaussian', U);
            case 't'
                [R,nu] = copulafit('t', U); %#ok<NASGU>
            otherwise
                theta = copulafit(fam, U); %#ok<NASGU>
        end
        V1 = U(:,1); V2 = nan(size(V1));
        for i=1:size(U,1)
            u1 = U(i,1); u2 = U(i,2);
            switch fam
                case 'Gaussian'
                    Cplus = copulacdf('Gaussian',[min(u1+epsFD,1) u2], R);
                    Cminus= copulacdf('Gaussian',[max(u1-epsFD,0) u2], R);
                case 't'
                    Cplus = copulacdf('t',[min(u1+epsFD,1) u2], R, nu);
                    Cminus= copulacdf('t',[max(u1-epsFD,0) u2], R, nu);
                otherwise
                    Cplus = copulacdf(fam,[min(u1+epsFD,1) u2], theta);
                    Cminus= copulacdf(fam,[max(u1-epsFD,0) u2], theta);
            end
            V2(i) = max(0,min(1,(Cplus - Cminus)/(2*epsFD)));
        end
        [~,p1] = kstest(V1,[sort(V1) ( (1:numel(V1))'/numel(V1) )]);
        [~,p2] = kstest(V2,[sort(V2) ( (1:numel(V2))'/numel(V2) )]);
        pR = 1 - abs(corr(V1,V2,'type','Spearman','rows','complete'));
        msg = sprintf('GOF (Rosenblatt 2D)\nKS p(V1~U): %.3f\nKS p(V2~U): %.3f\nIndependence proxy p: %.3f (higher ~ better)', p1,p2,pR);
        showAlert(src, msg, 'GOF results');
    end

%% ===================== Rolling Kendall tau =====================
    function rollingTauDialog(src)
        cols = getSelectedColumns(); if numel(cols)~=2, showAlert(src,'Select exactly 2 columns.','Rolling'); return; end
        dlg = uifigure('Name','Rolling tau settings','Position',[260 260 340 160]);
        uilabel(dlg,'Text','Window W:','Position',[20 90 100 22],'HorizontalAlignment','right');
        edW = uieditfield(dlg,'numeric','Position',[130 90 160 22],'Value',60,'Limits',[10 Inf],'RoundFractionalValues',true);
        uibutton(dlg,'push','Text','Run','Position',[110 40 120 26], 'ButtonPushedFcn', @(~,~) runRollingTau(dlg, round(edW.Value)));
    end

    function runRollingTau(dlg, W)
        try
            cols = getSelectedColumns(); X = userData.Table{:,cols}; marg = ddlMarg.Value; U = makeUfromMarginals(X, marg); n=size(U,1);
            tau = nan(n,1); for t=W:n, tau(t)=corr(U(t-W+1:t,1),U(t-W+1:t,2),'type','Kendall','rows','complete'); end
            figure('Name',sprintf('Rolling tau (W=%d)',W),'Color','w'); plot(tau,'-','LineWidth',1.3); grid on; ylim([-1 1]); yline(0,'k:');
            xlabel('t'); ylabel('Kendall''s tau'); title(sprintf('Rolling tau for %s-%s (W=%d)', cols{1}, cols{2}, W));
            try, close(dlg); end
        catch ME
            showAlert(dlg, ME.message, 'Rolling Error');
        end
    end

%% ===================== Forecast (GAS) =====================
    function forecastDialog(src)
        cols = getSelectedColumns(); if numel(cols)~=2, showAlert(src,'Select exactly 2 columns for forecasting (e.g., X1,X2).','Forecast'); return; end
        dlg = uifigure('Name','Forecast tau(t+1) settings','Position',[240 240 380 240]);
        uilabel(dlg,'Text','Family:','Position',[20 180 120 22],'HorizontalAlignment','right');
        ddFam = uidropdown(dlg,'Position',[150 180 200 22],'Items',{'Gaussian','t'},'Value','Gaussian');
        uilabel(dlg,'Text','W (rolling tau):','Position',[20 140 120 22],'HorizontalAlignment','right');
        edW = uieditfield(dlg,'numeric','Position',[150 140 200 22],'Value',60,'Limits',[10 Inf],'RoundFractionalValues',true);
        uilabel(dlg,'Text','\alpha (1-coverage):','Position',[20 100 120 22],'HorizontalAlignment','right');
        edA = uieditfield(dlg,'numeric','Position',[150 100 200 22],'Value',max(0.01,min(0.30,edtAlpha.Value)),'Limits',[0.01 0.3]);
        uilabel(dlg,'Text','Calibration split:','Position',[20 60 120 22],'HorizontalAlignment','right');
        edC = uieditfield(dlg,'numeric','Position',[150 60 200 22],'Value',0.20,'Limits',[0.05 0.9]);
        uibutton(dlg,'push','Text','Run Forecast','Position',[120 20 140 26], ...
            'ButtonPushedFcn', @(~,~)runForecast(dlg, ddFam.Value, round(edW.Value), edA.Value, edC.Value));
    end

    function runForecast(dlg, fam, W, alpha, calFrac)
        try
            cols = getSelectedColumns(); if numel(cols)~=2, showAlert(dlg,'Need exactly 2 columns.','Forecast'); return; end
            X = userData.Table{:,cols}; marg = ddlMarg.Value; out = forecastTauGAS_core(X, marg, W, fam, alpha, calFrac);
            n = size(X,1); tt = (1:n)'; m = ~isnan(out.tau_pred) & ~isnan(out.tau_real);
            fig = figure('Name',sprintf('tau Forecast (%s, W=%d, \alpha=%.2f)', fam, W, alpha),'Color','w'); %#ok<NASGU>
            plot(tt(m), out.tau_real(m), '-', 'LineWidth',1.2); hold on;
            plot(tt(m), out.tau_pred(m), '--', 'LineWidth',1.6);
            plot(tt(m), out.loA(m), '-.', 'LineWidth',1.0); plot(tt(m), out.hiA(m), '-.', 'LineWidth',1.0);
            plot(tt(m), out.lo90(m), ':', 'LineWidth',1.0); plot(tt(m), out.hi90(m), ':', 'LineWidth',1.0);
            plot(tt(m), out.lo95(m), ':', 'LineWidth',0.9); plot(tt(m), out.hi95(m), ':', 'LineWidth',0.9);
            yline(0,'k:'); grid on; xlabel('t'); ylabel('Kendall''s tau');
            legend('tau realized','tau^ (GAS)','(1-\alpha) lo','(1-\alpha) hi','90% lo','90% hi','95% lo','95% hi','Location','best');
            title(sprintf('%s + GAS | W=%d | \alpha=%.2f | cal=%.2f', fam, W, alpha, calFrac)); hold off;
            L = min(50, sum(m)); idx = find(m, L, 'last');
            fig2 = uifigure('Name','Forecast table','Position',[300 300 720 300]); %#ok<NASGU>
            uitable(fig2,'Data', [num2cell(idx(:)), ...
                                  num2cell(out.tau_real(idx)), num2cell(out.tau_pred(idx)), ...
                                  num2cell(out.loA(idx)), num2cell(out.hiA(idx)), ...
                                  num2cell(out.lo90(idx)), num2cell(out.hi90(idx)), ...
                                  num2cell(out.lo95(idx)), num2cell(out.hi95(idx))], ...
                           'ColumnName', {'t','tau_real','tau_pred','loA','hiA','lo90','hi90','lo95','hi95'}, ...
                           'Position',[20 20 680 250]);
            try, close(dlg); end
        catch ME
            showAlert(dlg, ME.message, 'Forecast Error');
        end
    end

%% ===================== Forecast core (nested) =====================
    function out = forecastTauGAS_core(X, marg, W, family, alpha, calFrac)
        n = size(X,1);
        if n < W+150
            error('Need at least ~%d observations (got %d).', W+150, n);
        end
        U = makeUfromMarginals(X, marg);
        Z = [norminv(U(:,1)), norminv(U(:,2))]; Z(~isfinite(Z))=0; Xex = abs(Z);
        tau_real = nan(n,1);
        for t=W:n, tau_real(t) = corr(U(t-W+1:t,1), U(t-W+1:t,2), 'type','Kendall','rows','complete'); end
        nu = [];
        if strcmpi(family,'t')
            try
                [~, nuHat] = copulafit('t', U(max(1,W):end,:)); nu = max(3, min(100, round(nuHat)));
            catch, nu = 8; end
        end
        idxStart = W+1; idxEnd = n-1; Ttot = idxEnd-idxStart+1; Tcal = max(50, round(calFrac*Ttot)); Ttrain = Ttot - Tcal; trIdx = idxStart:(idxStart+Ttrain-1);
        muX = mean(Xex(trIdx,:),1,'omitnan'); sX  = std(Xex(trIdx,:),0,1,'omitnan'); sX(sX==0)=1; XexS = (Xex - muX)./sX;
        tau0 = tau_real(idxStart-1);
        if isnan(tau0), tau0 = corr(U(max(1,idxStart-W):idxStart-1,1), U(max(1,idxStart-W):idxStart-1,2), 'type','Kendall','rows','complete'); end
        tau0 = clamp(tau0,-0.999,0.999); rho0 = sin(pi*tau0/2); f0 = atanh(clamp(rho0,-0.9999,0.9999));
        theta0 = [0, log(0.1), log(0.9/0.1), 0, 0];
        obj = @(th) gas_negloglik(th, U, XexS, trIdx, f0, family, nu);
        opts = optimset('Display','off','MaxFunEvals',5e4,'MaxIter',5e4);
        thetaHat = fminsearch(obj, theta0, opts); %#ok<NASGU>
        [~, outTrain] = gas_negloglik(thetaHat, U, XexS, trIdx, f0, family, nu); params = outTrain.params;
        [tau_pred, ~] = predict_tau(params, U, XexS, tau_real, idxStart, idxEnd, family, nu);
        mask = ~isnan(tau_pred) & ~isnan(tau_real); idxAll = find(mask & ( (1:n)'>=idxStart & (1:n)'<=idxEnd ));
        Tcal = max(50, round(calFrac*numel(idxAll))); calIdx = idxAll(end-Tcal+1:end);
        resid = abs(tau_real(calIdx) - tau_pred(calIdx)); resid = resid(isfinite(resid));
        if isempty(resid), qA=0.08; q90=0.10; q95=0.14; else, qA=quantile(resid,1-alpha); q90=quantile(resid,0.90); q95=quantile(resid,0.95); end
        loA  = clamp(tau_pred - qA,  -1, 1);  hiA  = clamp(tau_pred + qA,  -1, 1);
        lo90 = clamp(tau_pred - q90, -1, 1);  hi90 = clamp(tau_pred + q90, -1, 1);
        lo95 = clamp(tau_pred - q95, -1, 1);  hi95 = clamp(tau_pred + q95, -1, 1);
        out = struct('tau_real',tau_real,'tau_pred',tau_pred,'loA',loA,'hiA',hiA,'lo90',lo90,'hi90',hi90,'lo95',lo95,'hi95',hi95,'params',params,'nu',nu,'trainIdx',trIdx);
    end

%% ===================== GAS internals =====================
    function [NLL, out] = gas_negloglik(thetaU, U, Xex, trIdx, f0, fam, nu)
        % thetaU = [w, log(a), logit(b), g1, g2]
        w  = thetaU(1); a  = exp(thetaU(2)); b  = 1/(1+exp(-thetaU(3))); g1 = thetaU(4); g2 = thetaU(5); gam = [g1; g2];
        epsFD = 1e-4; n = size(U,1); f = nan(n,1); rho = nan(n,1); NLL = 0; f(trIdx(1)-1) = f0; rho(trIdx(1)-1) = tanh(f0);
        for t = trIdx(1):trIdx(end)
            xlag  = (t-1>=1) * Xex(max(t-1,1),:).'; fpred = w + b*f(t-1) + gam.'*xlag; rho_t = tanh(fpred);
            if ~isfinite(rho_t), rho_t = 0; end
            rho_t = max(min(rho_t, 0.9999), -0.9999);
            ll_t  = logCopulaPdf(U(t,:), rho_t, fam, nu); dldr  = dlogc_drho_fd(U(t,:), rho_t, fam, nu, epsFD); s_t = dldr * sech2(fpred);
            f(t)  = fpred + a * s_t; rho(t)= tanh(f(t)); NLL   = NLL - ll_t;
        end
        out = struct('params',struct('w',w,'a',a,'b',b,'g',gam),'f',f,'rho',rho);
    end

    function [tau_pred, fSeries] = predict_tau(params, U, Xex, tau_real, idxStart, idxEnd, fam, nu)
        n = size(U,1); w=params.w; a=params.a; b=params.b; gam=params.g; f = nan(n,1);
        tau0 = tau_real(idxStart-1);
        if isnan(tau0)
            tau0 = corr(U(max(1,idxStart-60):idxStart-1,1), U(max(1,idxStart-60):idxStart-1,2),'type','Kendall','rows','complete');
        end
        tau0 = clamp(tau0,-0.999,0.999); rho0 = sin(pi*tau0/2); f0 = atanh(clamp(rho0,-0.9999,0.9999)); f(idxStart-1)=f0;
        for t = idxStart:idxEnd
            xlag  = (t-1>=1) * Xex(max(t-1,1),:).'; f_pred= w + b*f(t-1) + gam.' * xlag; rho_t = tanh(f_pred);
            dldr  = dlogc_drho_fd(U(t,:), rho_t, fam, nu, 1e-4); s_t   = dldr * sech2(f_pred);
            f(t)  = f_pred + a * s_t;
        end
        tau_pred = nan(n,1);
        for t = idxStart:idxEnd
            xnow = Xex(t,:).'; f_for = w + b * f(t) + gam.' * xnow; rho_for = tanh(f_for);
            tau_pred(t+1) = (2/pi)*asin(clamp(rho_for,-0.9999,0.9999));
        end
        fSeries=f;
    end

    function val = logCopulaPdf(urow, rho, fam, nu)
        u = max(min(urow(:)', 1-1e-10), 1e-10);
        if ~isfinite(rho), rho = 0; end
        rho = max(min(rho, 0.9999), -0.9999);
        R = [1 rho; rho 1];
        if strcmpi(fam,'t')
            if isempty(nu), nu = 8; end
            c = copulapdf('t', u, R, nu);
        else
            c = copulapdf('Gaussian', u, R);
        end
        val = log(max(c, realmin));
    end

    function d = dlogc_drho_fd(urow, rho, fam, nu, epsFD)
        r1 = clamp(rho+epsFD, -0.9999, 0.9999); r2 = clamp(rho-epsFD, -0.9999, 0.9999);
        l1 = logCopulaPdf(urow, r1, fam, nu);  l2 = logCopulaPdf(urow, r2, fam, nu); d  = (l1 - l2) / (2*epsFD);
    end

    function y = sech2(x), y = 1./cosh(x).^2; end
    function z = clamp(x,a,b)
        if ~isfinite(x), z = 0; else, z = max(a, min(b, x)); end
    end

%% ===================== Advanced: K-fold CV =====================
    function kfoldCV(src)
        assertFitted(src);
        cols = getSelectedColumns(); if numel(cols)<2, showAlert(src,'Select ≥2','CV'); return; end
        X = userData.Table{:,cols};
        U = makeUfromMarginals(X, ddlMarg.Value);
        K = max(2, round(edtKfold.Value)); fams = {'Gaussian','t'};
        if size(U,2)==2, fams = [fams, {'Clayton','Frank','Gumbel'}]; end
        rng('default'); idx = crossvalind('Kfold', size(U,1), K);
        meanLL = nan(1,numel(fams));
        for fi=1:numel(fams)
            fam = fams{fi}; ll = zeros(K,1);
            for k=1:K
                tr = idx~=k; te = idx==k; Utr = U(tr,:); Ute = U(te,:);
                try
                    switch fam
                        case 'Gaussian'
                            R = copulafit('Gaussian', Utr);
                            ll(k) = sum(log(max(copulapdf('Gaussian', Ute, R), realmin)));
                        case 't'
                            [R,nu] = copulafit('t', Utr);
                            ll(k) = sum(log(max(copulapdf('t', Ute, R, max(2,round(nu))), realmin)));
                        otherwise
                            theta = copulafit(fam, Utr);
                            ll(k) = sum(log(max(copulapdf(fam, Ute, theta), realmin)));
                    end
                catch
                    ll(k) = -Inf;
                end
            end
            meanLL(fi) = mean(ll);
        end
        figure('Name','K-fold CV Log-Likelihood'); bar(meanLL);
        set(gca,'XTickLabel',fams,'XTick',1:numel(fams)); ylabel('Mean held-out loglik'); grid on;
        title(sprintf('K-fold (K=%d) CV log-likelihood', K));
    end

%% ===================== Advanced: Bootstrap =====================
    function bootstrapDialog(src)
        assertFitted(src);
        B = max(50, round(edtBootB.Value));
        choice = questdlg(sprintf('Run parametric bootstrap with B=%d?',B), 'Bootstrap','Run','Cancel','Run');
        if ~strcmpi(choice,'Run'), return; end
        bootstrapRun(src, B);
    end

    function bootstrapRun(src, B)
        assertFitted(src); fam = userData.family; U = userData.U; d=size(U,2);
        switch fam
            case {'Gaussian','t'}
                Rhat = userData.R; nu = userData.nu; if isempty(Rhat), showAlert(src,'No R','Bootstrap'); return; end
                Rhats = zeros([size(Rhat) B]); nuv = nan(B,1);
                for b=1:B
                    try
                        if strcmpi(fam,'Gaussian')
                            Ub = copularnd('Gaussian', Rhat, size(U,1));
                            Rb = copulafit('Gaussian', Ub);
                            Rhats(:,:,b) = Rb;
                        else
                            Ub = copularnd('t', Rhat, size(U,1), max(2,round(nu)));
                            [Rb, nub] = copulafit('t', Ub);
                            Rhats(:,:,b) = Rb; nuv(b)=nub;
                        end
                    catch
                        Rhats(:,:,b) = nan(size(Rhat)); nuv(b)=nan;
                    end
                end
                pairs = nchoosek(1:d,2); out = cell(size(pairs,1),5);
                for i=1:size(pairs,1)
                    a=pairs(i,1); b=pairs(i,2);
                    v = squeeze(Rhats(a,b,:)); v=v(isfinite(v));
                    if isempty(v), lo=NaN; hi=NaN; med=NaN; else
                        pr = prctile(v,[2.5 50 97.5]); lo=pr(1); med=pr(2); hi=pr(3);
                    end
                    out(i,:)={sprintf('R(%d,%d)',a,b), userData.R(a,b), lo, med, hi};
                end
                fig = uifigure('Name','Bootstrap CI for R','Position',[300 300 640 280]);
                uitable(fig,'Data',out,'ColumnName',{'Entry','Rhat','lo2.5%','median','hi97.5%'},'Position',[20 20 600 240]);
                if strcmpi(fam,'t')
                    nv = nuv(isfinite(nuv)); if ~isempty(nv)
                        pr = prctile(nv,[2.5 50 97.5]);
                        showAlert(src,sprintf('nu: %.2f (2.5%%=%.2f, 97.5%%=%.2f)',userData.nu,pr(1),pr(3)),'Bootstrap nu');
                    end
                end
            otherwise
                theta = userData.theta; if isempty(theta), showAlert(src,'No theta','Bootstrap'); return; end
                TH = nan(B,1);
                for b=1:B
                    try
                        Ub = copularnd(fam, theta, size(U,1));
                        thb = copulafit(fam, Ub); TH(b)=thb;
                    catch
                        TH(b)=nan;
                    end
                end
                v = TH(isfinite(TH)); if isempty(v), pr=[NaN NaN NaN]; else, pr = prctile(v,[2.5 50 97.5]); end
                showAlert(src, sprintf('theta: %.4f (2.5%%=%.4f, 97.5%%=%.4f)', theta, pr(1), pr(3)), 'Bootstrap theta');
        end
    end

%% ===================== Advanced: 3D scatter =====================
    function Uscatter3D(src)
        if isempty(userData.Table), showAlert(src,'Load data first','3D'); return; end
        cols = getSelectedColumns();
        if numel(cols) >= 3
            X = userData.Table{:,cols(1:3)};
            figure('Name','3D Scatter'); scatter3(X(:,1),X(:,2),X(:,3),20,X(:,3),'filled'); grid on;
            xlabel(cols{1}); ylabel(cols{2}); zlabel(cols{3}); title('3D Scatter (raw)'); colorbar; return;
        elseif numel(cols) == 2
            try
                figp = ancestor(src,'figure'); if isempty(figp) || ~isvalid(figp), figp = f; end
                choice = uiconfirm(figp,'Use Z = copula density c(u) (needs fit), or Z = empirical copula C_n(u)?', ...
                                   'Choose Z', 'Options',{'c(u) (fitted)','C_n(u)'}, ...
                                   'DefaultOption',1,'CancelOption',2);
            catch
                choice = 'C_n(u)';
            end
            X = userData.Table{:,cols}; U = makeUfromMarginals(X, ddlMarg.Value);
            U1 = U(:,1); U2 = U(:,2); Z = [];
            if strcmp(choice,'c(u) (fitted)')
                if isempty(userData.family)
                    showAlert(src,'Fit a copula first, then try again (for c(u)).','3D'); return;
                end
                fam = userData.family;
                try
                    switch fam
                        case 'Gaussian'
                            R = copulafit('Gaussian',[U1 U2]); Z = copulapdf('Gaussian',[U1 U2],R);
                        case 't'
                            [R,nu] = copulafit('t',[U1 U2]); Z = copulapdf('t',[U1 U2],R,max(2,round(nu)));
                        otherwise
                            th = copulafit(fam,[U1 U2]); Z = copulapdf(fam,[U1 U2],th);
                    end
                catch
                    showAlert(src,'Could not compute copula density. Falling back to C_n(u).','3D');
                    choice = 'C_n(u)';
                end
            end
            if isempty(Z) || strcmp(choice,'C_n(u)')
                n = numel(U1); Z = zeros(n,1);
                for i=1:n, Z(i) = mean( (U1<=U1(i)) & (U2<=U2(i)) ); end
            end
            figure('Name','U-Scatter 3D (derived Z)'); scatter3(U1,U2,Z,18,Z,'filled'); grid on; xlim([0 1]); ylim([0 1]);
            xlabel('U1'); ylabel('U2'); if strcmp(choice,'C_n(u)'), zlabel('Z = C_n(u)'); else, zlabel('Z = c(u)'); end
            title(sprintf('3D on U for %s-%s', cols{1}, cols{2})); colorbar; return;
        else
            showAlert(src,'Select 2 or 3 columns.','3D'); return;
        end
    end

%% ===================== Advanced: Conditional slice =====================
    function conditionalSliceDialog(src)
        if isempty(userData.Table), showAlert(src,'Load data first','Conditional'); return; end
        cols = getSelectedColumns(); if numel(cols)~=2, showAlert(src,'Select exactly 2 columns','Conditional'); return; end
        dlg = uifigure('Name','Conditional slice U2|U1=u','Position',[260 260 380 180]);
        uilabel(dlg,'Text','u (conditioning value for U1):','Position',[20 100 220 22]);
        edU = uieditfield(dlg,'numeric','Position',[250 100 100 22],'Value',0.5,'Limits',[0.01 0.99]);
        uibutton(dlg,'push','Text','Plot','Position',[140 50 100 26], 'ButtonPushedFcn', @(~,~) conditionalSlice(edU.Value, cols, dlg));
    end

    function conditionalSlice(uval, cols, dlg)
        try
            X = userData.Table{:,cols};
            U = makeUfromMarginals(X, ddlMarg.Value);
            fam = userData.family;
            switch fam
                case 'Gaussian'
                    R = copulafit('Gaussian', U);
                case 't'
                    [R,nu] = copulafit('t', U);
                otherwise
                    R = []; nu = []; % not needed for Archimedean
            end
            v = linspace(0.001,0.999,400)'; u = max(0.001,min(0.999,uval));
            switch fam
                case 'Gaussian'
                    rho = R(1,2);
                    z1 = norminv(u); mu = rho*z1; s2 = max(1e-8,1-rho^2);
                    z2 = norminv(v);
                    fz = normpdf(z2, mu, sqrt(s2)); fu = normpdf(z2);
                    cond = fz ./ fu;
                case 't'
                    rho = R(1,2); nu0 = max(2,round(nu));
                    z1 = tinv(u, nu0);
                    mu = rho*z1; s2 = (nu0+z1.^2).*(1-rho^2)/(nu0+1);
                    z2 = tinv(v, nu0+1);
                    fz = tpdf((z2-mu)./sqrt(s2), nu0+1) ./ sqrt(s2);
                    fu = tpdf(z2, nu0+1); cond = fz ./ fu;
                otherwise
                    th = copulafit(fam, U);
                    Cline = copulapdf(fam, [repmat(u,numel(v),1) v], th);
                    cond = Cline / trapz(v, Cline);
            end
            figure('Name','Conditional slice U2|U1=u'); plot(v, cond, 'LineWidth',1.4); grid on; xlim([0 1]);
            xlabel(sprintf('v = U_2 | U_1 = %.2f',u)); ylabel('density'); title(sprintf('Conditional density slice (%s)', fam));
            try, close(dlg); end
        catch ME
            showAlert(dlg, ME.message, 'Conditional Error');
        end
    end

%% ===================== Simulation =====================
    function simulateFromCopula(src)
        if isempty(userData.U) || (isempty(userData.R) && isempty(userData.theta))
            showAlert(src,'Fit a copula first','Error'); return;
        end
        cols = getSelectedColumns(); d = numel(cols);
        n = max(10, round(edtNSim.Value)); marg = ddlMarg.Value;
        switch lower(userData.family)
            case 't'
                Usim = copularnd('t', userData.R, n, max(2,round(userData.nu)));
            case 'gaussian'
                Usim = copularnd('Gaussian', userData.R, n);
            case 'clayton'
                if d~=2, showAlert(src,'Archimedean simulation requires exactly 2 columns.','Error'); return; end
                Usim = copularnd('Clayton', userData.theta, n);
            case 'frank'
                if d~=2, showAlert(src,'Archimedean simulation requires exactly 2 columns.','Error'); return; end
                Usim = copularnd('Frank', userData.theta, n);
            case 'gumbel'
                if d~=2, showAlert(src,'Archimedean simulation requires exactly 2 columns.','Error'); return; end
                Usim = copularnd('Gumbel', userData.theta, n);
            otherwise
                showAlert(src,'Unknown family','Error'); return;
        end
        X_sim = zeros(n,d);
        for i = 1:d
            x_real = userData.Table{:,cols{i}};
            if strcmpi(marg,'Ranks (empirical)')
                X_sim(:,i) = empirical_icdf(x_real, Usim(:,i));
            else
                try
                    pd = fitdist(x_real, lower(marg));
                    X_sim(:,i) = icdf(pd, Usim(:,i));
                catch
                    X_sim(:,i) = empirical_icdf(x_real, Usim(:,i));
                end
            end
        end
        figure('Name','Simulated Data'); plotmatrix(X_sim); sgtitle(sprintf('Simulated Data from %s copula', userData.family));
    end

    function plotJointPDF(src)
        cols = getSelectedColumns();
        if numel(cols) ~= 2, showAlert(src, 'Select exactly 2 columns for joint PDF plot.', 'Invalid Selection'); return; end
        if isempty(userData.family), showAlert(src,'Fit a copula first.','Error'); return; end
        X1 = userData.Table{:,cols{1}}; X2 = userData.Table{:,cols{2}}; marg = ddlMarg.Value;
        U = makeUfromMarginals([X1 X2], marg); U1=U(:,1); U2=U(:,2);
        famKey = userData.family; R2=[]; nu=[]; th=[];
        switch famKey
            case 't'
                [R2,nu]=copulafit('t',[U1 U2]);
            case 'Gaussian'
                R2=copulafit('Gaussian',[U1 U2]);
            case {'Clayton','Frank','Gumbel'}
                th = copulafit(famKey,[U1 U2]);
            otherwise
                showAlert(src,'Unsupported for PDF plot','Error'); return;
        end
        x1 = linspace(min(X1), max(X1), 80); x2 = linspace(min(X2), max(X2), 80);
        [Xgrid, Ygrid] = meshgrid(x1, x2);
        if strcmpi(marg,'Ranks (empirical)')
            U1g = empirical_cdf(X1, Xgrid); U2g = empirical_cdf(X2, Ygrid);
            f1  = empirical_pdf(X1, Xgrid); f2  = empirical_pdf(X2, Ygrid);
        else
            pd1 = fitdist(X1, lower(marg)); pd2 = fitdist(X2, lower(marg));
            U1g = cdf(pd1, Xgrid); U2g = cdf(pd2, Ygrid);
            f1  = pdf(pd1, Xgrid); f2  = pdf(pd2, Ygrid);
        end
        switch famKey
            case 't'
                c = copulapdf('t',[U1g(:) U2g(:)], R2, max(2,round(nu)));
            case 'Gaussian'
                c = copulapdf('Gaussian',[U1g(:) U2g(:)], R2);
            otherwise
                c = copulapdf(famKey,[U1g(:) U2g(:)], th);
        end
        c = reshape(c, size(Xgrid)); jointPDF = c .* f1 .* f2;
        figure('Name','Joint PDF'); surf(x1, x2, jointPDF, 'EdgeColor','none');
        xlabel(cols{1}); ylabel(cols{2}); zlabel('Joint PDF');
        title(sprintf('Joint PDF: %s + %s', userData.family, marg)); view(135,30); colorbar;
    end

   function autoSelectCopula(src)
        cols = getSelectedColumns();
        if isempty(cols) || numel(cols) < 2
            showAlert(src,'Select at least 2 columns','Error');
            return;
        end
        X = userData.Table{:,cols};
        marg = ddlMarg.Value;
        U = makeUfromMarginals(X, marg);
        families = {'Gaussian','t'};
        if size(U,2)==2
            families = [families, {'Clayton','Frank','Gumbel'}];
        end
        results = {};
        for k = 1:numel(families)
            try
                fam = families{k};
                if strcmpi(fam,'t')
                    [R, nu] = copulafit('t', U);
                    logL = sum(log(max(copulapdf('t', U, R, nu), realmin)));
                    kparams = numel(R(triu(true(size(R)),1))) + 1;
                elseif strcmpi(fam,'gaussian')
                    R = copulafit('Gaussian', U);
                    logL = sum(log(max(copulapdf('Gaussian', U, R), realmin)));
                    kparams = numel(R(triu(true(size(R)),1)));
                else
                    theta = copulafit(fam, U); %#ok<NASGU>
                    logL = sum(log(max(copulapdf(fam, U, theta), realmin)));
                    kparams = 1;
                end
                nObs = size(U,1);
                AIC = -2*logL + 2*kparams;
                BIC = -2*logL + kparams*log(nObs);
                results(end+1,:) = {fam, logL, AIC, BIC}; %#ok<AGROW>
            catch
                results(end+1,:) = {families{k}, NaN, NaN, NaN}; %#ok<AGROW>
            end
        end
        T = cell2table(results, 'VariableNames', {'Copula','LogL','AIC','BIC'});
        [~, iAIC] = min(T.AIC);
        [~, iBIC] = min(T.BIC);
        msg = sprintf('Model selection completed:\nBest AIC: %s\nBest BIC: %s', ...
                      T.Copula{iAIC}, T.Copula{iBIC});
        showAlert(src, msg, 'Copula Selection Results');
        figResults = uifigure('Name','Model Selection Table','Position',[320 320 520 220]);
        data = [string(T.Copula), string(round(T.LogL,3)), string(round(T.AIC,3)), string(round(T.BIC,3))];
        uitable(figResults, 'Data', data, 'ColumnName', {'Copula','LogL','AIC','BIC'}, 'Position', [20 20 480 170]);
    end

%% ===================== Export / Session =====================
    function exportResults(src)
        if isempty(userData.U), showAlert(src,'Nothing to export (fit first).','Export'); return; end
        [file,path] = uiputfile('copula_results.mat','Save Results'); if isequal(file,0), return; end
        S = struct('U',userData.U,'R',userData.R,'nu',userData.nu,'theta',userData.theta,'family',userData.family,'marg',userData.marg,'Cluster',userData.Cluster);
        save(fullfile(path,file),'-struct','S'); showAlert(src,'Saved results MAT file.','Export');
    end

    function saveSession(src)
        if isempty(userData.Table), showAlert(src,'No session to save.','Save'); return; end
        [file,path] = uiputfile('copula_session_allinone.mat','Save Session'); if isequal(file,0), return; end
        S = struct('Table',userData.Table,'U',userData.U,'R',userData.R,'nu',userData.nu,'theta',userData.theta,'family',userData.family,'marg',userData.marg, 'Cluster',userData.Cluster);
        save(fullfile(path,file),'-struct','S'); showAlert(src,'Session saved.','Save');
    end

    function loadSession(src)
        [file,path] = uigetfile('*.mat','Load Session'); if isequal(file,0), return; end
        S = load(fullfile(path,file));
        if isfield(S,'Table')
            userData.Table = S.Table; names = userData.Table.Properties.VariableNames;
            lstCols.Items = names; if numel(names)>=2, lstCols.Value = names(1:2); else, lstCols.Value=names(1); end
        end
        fields = {'U','R','nu','theta','family','marg','Cluster'};
        for k=1:numel(fields), if isfield(S,fields{k}), userData.(fields{k}) = S.(fields{k}); end, end
        txtOut.Value = 'Session loaded.';
        if ~isempty(userData.R)
            imagesc(ax,userData.R); colorbar(ax); title(ax,'Copula Correlation (R)');
        end
    end

    function exportAxesPNG(src)
        try
            [file,path] = uiputfile('figure.png','Export current axes'); if isequal(file,0), return; end
            exportgraphics(ax, fullfile(path,file)); showAlert(src,'PNG saved.','Export');
        catch ME
            showAlert(src,ME.message,'Export Error');
        end
    end

    function quickReport(src)
        try
            ts = datestr(now,'yyyymmdd_HHMMSS');
            d = uigetdir(pwd, 'Select folder to create report'); if isequal(d,0), return; end
            rep = fullfile(d, ['report_' ts]); if ~exist(rep,'dir'), mkdir(rep); end
            try, exportgraphics(ax, fullfile(rep,'main_axes.png')); end
            fid = fopen(fullfile(rep,'summary.txt'),'w');
            fprintf(fid,'Copula Toolbox Quick Report\n');
            fprintf(fid,'-----------------------------\n');
            fprintf(fid,'Family: %s\nMarginals: %s\n', userData.family, userData.marg);
            if ~isempty(userData.R)
                fprintf(fid,'R matrix:\n'); fprintf(fid,'%s\n', evalc('disp(userData.R)'));
            end
            if ~isempty(userData.nu), fprintf(fid,'nu: %.4f\n', userData.nu); end
            if ~isempty(userData.theta), fprintf(fid,'theta: %.4f\n', userData.theta); end
            if isfield(userData,'Cluster') && ~isempty(userData.Cluster) && isfield(userData.Cluster,'results') && ~isempty(userData.Cluster.results)
                CC = userData.Cluster.results;
                fprintf(fid,'\n[Cluster] Alg=%s | K=%d | UseU=%d | Marg=%s\n', CC.alg, CC.K, userData.Cluster.lastUseU, CC.marginalName);
                for kk=1:CC.K
                    fprintf(fid,'  - Cluster %d: family=%s | size=%d | BIC=%.1f\n', kk, CC.clusters(kk).fit.family, CC.clusters(kk).size, CC.clusters(kk).criteria.BIC);
                end
            end
            fclose(fid);
            showAlert(src, sprintf('Report created at\n%s', rep), 'Report');
        catch ME
            showAlert(src, ME.message, 'Report Error');
        end
    end

%% ===================== Utilities =====================
    function A = nearestPD(A)
        [V,D] = eig((A+A')/2); D = diag(D); D(D<1e-10) = 1e-10; A = V*diag(D)*V'; A = (A+A')/2;
        s = sqrt(diag(A)); A = A ./ (s*s'); A(1:size(A,1)+1:end) = 1;
    end

    function q = empirical_icdf(x, u)
        x = x(:); u = u(:); [xs,~] = sort(x);
        ranks = (1:numel(x))'/(numel(x)+1);
        q = interp1(ranks, xs, u, 'linear','extrap');
    end

    function Fu = empirical_cdf(x, grid)
        x = x(:); xs = sort(x);
        Fu = arrayfun(@(z) mean(xs<=z, 'omitnan'), grid);
    end

    function fx = empirical_pdf(x, grid)
        x = x(:); if numel(x)<5, fx = zeros(size(grid)); return; end
        bw = 1.06*std(x,'omitnan')*numel(x)^(-1/5);
        if bw<=0 || ~isfinite(bw), bw = max(eps, iqr(x)/1.34* numel(x)^(-1/5)); end
        fx = zeros(size(grid));
        for i=1:numel(x)
            fx = fx + (1/(sqrt(2*pi)*bw)) * exp(-0.5*((grid - x(i))/bw).^2);
        end
        fx = fx / numel(x);
    end

    function showAlert(src,msg,title)
        if nargin<3 || isempty(title), title = 'Info'; end
        try
            fig = [];
            if ~isempty(src)
                try, if isvalid(src), fig = ancestor(src,'figure'); end, catch, fig = []; end
            end
            if isempty(fig) || ~isvalid(fig)
                try, if exist('f','var') && isvalid(f), fig = f; end, catch, fig = []; end
            end
            if isempty(fig) || ~isvalid(fig), fig = uifigure('Name',title); end
            uialert(fig, msg, title);
        catch
            try
                warndlg(msg, title);
            catch
                fprintf(2,'%s: %s\n', char(title), char(msg));
            end
        end
    end

%% ===================== Cluster→Copula Pipeline =====================
    function runClusterAndFit(src)
        if isempty(userData.Table)
            showAlert(src,'Load data first','Cluster'); return;
        end
        cols = getSelectedColumns();
        if numel(cols) < 2
            showAlert(src,'Select ≥2 columns','Cluster'); return;
        end
        X = userData.Table{:,cols};
        marg = ddlMarg.Value;

        alg = ddlClAlg.Value;
        Ktxt = strtrim(edtClK.Value);
        if isempty(Ktxt), Kset = []; else, Kset = str2double(Ktxt); if isnan(Kset) || Kset<2, Kset=[]; end, end
        Kmax = max(2, round(edtClKmax.Value));
        useU = chkClUseU.Value;

        epsv = edtClEps.Value; if ~isfinite(epsv) || epsv<=0, epsv = []; end
        minv = edtClMin.Value; if ~isfinite(minv) || minv<=0, minv = []; end
        fitMix = chkMix2D.Value; mixMax = round(edtMixK.Value);

        try
            CC = cluster_copula_module_inline(X, ...
                    'Marginals',marg, 'Algorithm',alg, 'K',Kset, 'Kmax',Kmax, ...
                    'UseUspace',useU, 'Epsilon',epsv, 'MinPts',minv, ...
                    'FitMixture2D',fitMix, 'MixMaxComp',mixMax, ...
                    'CopulaFamilies', {'Gaussian','t','Clayton','Frank','Gumbel'}, ...
                    'Verbose', true);
            % --- Backward-compatibility aliases (ώστε τα dialogs να "βλέπουν" σωστά) ---
CC.labels    = CC.clust;
CC.algorithm = CC.alg;
CC.useU      = useU;

        catch ME
            showAlert(src, ME.message, 'Cluster Error');
            return;
        end

        userData.Cluster.results = CC;
        userData.Cluster.lastAlg = alg; userData.Cluster.lastK = CC.K; userData.Cluster.lastUseU = useU;

                s = sprintf('Cluster→Copula done.\nAlg=%s | K=%d | UseU=%d | Marg=%s\n', ...
            alg, CC.K, useU, marg);
        for k = 1:CC.K
            cfit = CC.clusters(k).fit;
            if isfield(cfit,'R') && ~isempty(cfit.R)
                famShow = [upper(cfit.family(1)) cfit.family(2:end)];
                s = [s, sprintf('  Cluster %d: %s  | BIC=%.1f\n', k, famShow, CC.clusters(k).criteria.BIC)];
            else
                if isfield(cfit,'theta') && ~isempty(cfit.theta)
                    s = [s, sprintf('  Cluster %d: %s(theta=%.3f) | BIC=%.1f\n', k, upper(cfit.family), cfit.theta, CC.clusters(k).criteria.BIC)];
                else
                    s = [s, sprintf('  Cluster %d: %s | BIC=%.1f\n', k, upper(cfit.family), CC.clusters(k).criteria.BIC)];
                end
            end
        end
txtOut.Value = char(s);
        showAlert(src,'Clustering + per-cluster copulas completed.','Cluster');
    end

    function showClusterPlots(src)
        if isempty(userData.Cluster) || ~isfield(userData.Cluster,'results') || isempty(userData.Cluster.results)
            showAlert(src,'Run Cluster & Fit first','Plots'); return;
        end
        CC = userData.Cluster.results;
        try
            if isfield(CC.plots,'clusterScatter') && ~isempty(CC.plots.clusterScatter) && isvalid(CC.plots.clusterScatter)
                figure(CC.plots.clusterScatter);
            end
        catch, end
        try
            if isfield(CC.plots,'clusterUContours') && ~isempty(CC.plots.clusterUContours) && isvalid(CC.plots.clusterUContours)
                figure(CC.plots.clusterUContours);
            end
        catch, end
        if size(CC.U,2)==2 && ~isempty(CC.mixture2D)
            Rsp = CC.mixture2D.Rsp; [~,lab] = max(Rsp,[],2);
            figure('Name','Mixture-of-Copulas responsibilities','Color','w');
            scatter(CC.U(:,1),CC.U(:,2),12,lab,'filled'); xlim([0 1]); ylim([0 1]); colorbar; grid on;
            title(sprintf('Mixture-of-Copulas (K=%d)', CC.mixture2D.K));
        end
    end

    function exportClusterResults(src)
        if isempty(userData.Cluster) || ~isfield(userData.Cluster,'results') || isempty(userData.Cluster.results)
            showAlert(src,'Nothing to export','Export Cluster'); return;
        end
        [file,path] = uiputfile('cluster_copula_results.mat','Save Cluster Results');
        if isequal(file,0), return; end
        CC = userData.Cluster.results;
        save(fullfile(path,file),'-struct','CC');
        showAlert(src,'Cluster results saved.','Export');
    end



%% ===================== Extra Analytics =====================
    function hsicTestDialog(src)
        cols = getSelectedColumns(); 
        if numel(cols)~=2, showAlert(src,'Select exactly 2 columns for HSIC.','HSIC'); return; end
        X = userData.Table{:,cols};
        U = makeUfromMarginals(X, ddlMarg.Value);
        [stat,p] = hsic_test(U(:,1), U(:,2));
        msg = sprintf('HSIC independence test\nstat=%.4f, p≈%.4f (small p => dependence)', stat, p);
        showAlert(src, msg, 'HSIC');
        figure('Name','HSIC Kernel Density'); 
        [k1,~]=ksdensity(U(:,1)); [k2,~]=ksdensity(U(:,2));
        plot(k1,'-'); hold on; plot(k2,'--'); legend('U1','U2'); grid on; title('Univariate densities (for context)');
    end

    function [stat,p] = hsic_test(x,y)
        x=x(:); y=y(:);
        n = numel(x);
        sigma_x = median(pdist(x)); if ~isfinite(sigma_x) || sigma_x<=0, sigma_x=0.1; end
        sigma_y = median(pdist(y)); if ~isfinite(sigma_y) || sigma_y<=0, sigma_y=0.1; end
        K = rbf_kernel(x,x,sigma_x); L = rbf_kernel(y,y,sigma_y);
        H = eye(n) - ones(n)/n;
        Kc = H*K*H; Lc = H*L*H;
        stat = trace(Kc*Lc)/(n-1)^2;
        % gamma approximation (Gretton et al.)
        mu = mean(Kc(:))*mean(Lc(:));
        varK = var(Kc(:)); varL = var(Lc(:)); v = varK*varL;
        if v<=0, v=eps; end
        alpha = mu^2 / v; beta = v / mu;
        p = 1 - gamcdf(stat, alpha, beta);
        if ~isfinite(p) || p<0, p=0; elseif p>1, p=1; end
    end

    function K = rbf_kernel(a,b,s)
        D = pdist2(a,b).^2;
        K = exp(-D/(2*s^2));
    end

    function covarDialog(src)
        cols = getSelectedColumns(); 
        if numel(cols)~=2, showAlert(src,'Select exactly 2 columns.','CoVaR'); return; end
        if isempty(userData.family), showAlert(src,'Fit a copula first.','CoVaR'); return; end
        X = userData.Table{:,cols}; 
        alpha = max(0.01,min(0.2, edtAlpha.Value));
        % compute CoVaR of X2 given X1 in distress (VaR_alpha)
        [VaR1, CoVaR2, CoES2] = copula_covar_coes(X(:,1), X(:,2), userData.marg, userData.family, userData.R, userData.nu, userData.theta, alpha);
        msg = sprintf('alpha=%.2f\nVaR_%s(X1)=%.4g\nCoVaR_%s(X2|X1 distress)=%.4g\nCoES_%s(X2|X1 distress)=%.4g', ...
                       alpha, num2str(alpha), VaR1, num2str(alpha), CoVaR2, num2str(alpha), CoES2);
        showAlert(src, msg, 'CoVaR/CoES');
    end

    function [VaR1, CoVaR2, CoES2] = copula_covar_coes(X1, X2, marg, fam, R, nu, theta, alpha)
        % map marginals to U
        U = makeUfromMarginals([X1 X2], marg); U1=U(:,1); U2=U(:,2);
        % VaR of X1
        if strcmpi(marg,'Ranks (empirical)')
            VaR1 = empirical_icdf(X1, alpha);
        else
            try, pd1 = fitdist(X1, lower(marg)); VaR1 = icdf(pd1, alpha); catch, VaR1 = quantile(X1, alpha); end
        end
        % conditional distribution of U2 | U1 = u_alpha
        ua = mean(U1 <= alpha);
        vgrid = linspace(0.001,0.999,200)';
        switch fam
            case 'Gaussian'
                rho = R(1,2); z1 = norminv(ua); mu = rho*z1; s2 = max(1e-8,1-rho^2);
                z2 = norminv(vgrid); fz = normpdf(z2, mu, sqrt(s2)); fu = normpdf(z2); cond = fz./fu;
            case 't'
                rho = R(1,2); nu0 = max(2,round(nu)); z1 = tinv(ua, nu0);
                mu = rho*z1; s2 = (nu0+z1.^2).*(1-rho^2)/(nu0+1);
                z2 = tinv(vgrid, nu0+1); fz = tpdf((z2-mu)./sqrt(s2), nu0+1) ./ sqrt(s2); fu = tpdf(z2, nu0+1); cond = fz./fu;
            otherwise
                th = theta; Cline = copulapdf(fam, [repmat(ua,numel(vgrid),1) vgrid], th); cond = Cline / trapz(vgrid, Cline);
        end
        Fcond = cumtrapz(vgrid, cond); Fcond = Fcond / max(Fcond);
        % CoVaR on U2 at alpha
        u2star = interp1(Fcond, vgrid, alpha, 'linear','extrap'); u2star = max(0,min(1,u2star));
        % map back to X2
        if strcmpi(marg,'Ranks (empirical)')
            CoVaR2 = empirical_icdf(X2, u2star);
            CoES2  = mean(X2(U2<=u2star));
        else
            try, pd2 = fitdist(X2, lower(marg)); CoVaR2 = icdf(pd2, u2star); CoES2 = mean(icdf(pd2, U2(U2<=u2star))); 
            catch, CoVaR2 = quantile(X2, u2star); CoES2 = mean(X2(X2<=CoVaR2)); end
        end
    end

    function clusterStabilityDialog(src)
        if isempty(userData.Cluster) || ~isfield(userData.Cluster,'results') || isempty(userData.Cluster.results)
            showAlert(src,'Run Cluster & Fit first','Stability'); return;
        end
        CC = userData.Cluster.results;
        B = max(50, round(edtBootB.Value));
        stab = cluster_consistency(CC.X, CC.clust, B);
        figure('Name','Cluster Stability','Color','w');
        bar(stab.per_cluster_jaccard); grid on; xlabel('Cluster'); ylabel('Avg Jaccard');
        title(sprintf('Bootstrap stability (B=%d)',B));
        msg = sprintf('Overall mean Jaccard = %.3f', stab.mean_jaccard);
        showAlert(src, msg, 'Stability');
    end

    function out = cluster_consistency(X, labels, B)
        K = max(labels); n = size(X,1);
        J = zeros(B, K);
        for b=1:B
            idx = randsample(n, n, true);
            Xb = X(idx,:); lab0 = labels(idx);
            labb = kmeans(Xb, K, 'Replicates',5,'MaxIter',300,'Display','off');
            for k=1:K
                set0 = (lab0==k); set1 = (labb==k);
                J(b,k) = jaccard_idx(set0, set1);
            end
        end
        out.per_cluster_jaccard = mean(J,1,'omitnan');
        out.mean_jaccard = mean(out.per_cluster_jaccard,'omitnan');
    end

    function j = jaccard_idx(a,b)
        a = a(:)>0; b=b(:)>0; j = sum(a&b) / max(1,sum(a|b));
    end

    function elbowSilhouette(src)
        cols = getSelectedColumns();
        if numel(cols)<2, showAlert(src,'Select ≥2 columns','Elbow'); return; end
        X = userData.Table{:,cols}; 
        useU = chkClUseU.Value;
        Z = X; if useU, Z = makeUfromMarginals(X, ddlMarg.Value); end
        Kmax = max(3, round(edtClKmax.Value)); 
        Kcand = 2:Kmax;
        SSE = zeros(size(Kcand)); SIL = zeros(size(Kcand));
        for ii=1:numel(Kcand)
            K = Kcand(ii);
            [lab, C, sumd] = kmeans(Z, K, 'Replicates',5,'MaxIter',500,'Display','off'); %#ok<ASGLU>
            SSE(ii) = sum(sumd);
            try
                SIL(ii) = mean(silhouette(Z, lab));
            catch
                SIL(ii) = NaN;
            end
        end
        figure('Name','Elbow & Silhouette','Color','w');
        yyaxis left; plot(Kcand, SSE,'-o'); ylabel('SSE');
        yyaxis right; plot(Kcand, SIL,'-s'); ylabel('Mean silhouette');
        xlabel('K'); grid on; title('Elbow & Silhouette diagnostics');
    end

%% ===================== Inline Cluster Module =====================
    function CC = cluster_copula_module_inline(X, varargin)
        % Ενσωματωμένη εκδοχή του module (ίδιο API)
        ip = inputParser; ip.FunctionName = 'cluster_copula_module_inline';
        addRequired(ip,'X', @(z)isnumeric(z) && ismatrix(z) && size(z,2)>=2);
        validMarg = {'Normal','Lognormal','Exponential','Gamma','tLocationScale','Ranks (empirical)'};
        addParameter(ip,'Marginals','Normal', @(s)ischar(s) && any(strcmpi(s,validMarg)));
        addParameter(ip,'Algorithm','kmeans', @(s)ischar(s) && any(strcmpi(s,{'kmeans','gmm','hier','dbscan'})));
        addParameter(ip,'K',[], @(z)isempty(z) || (isscalar(z) && z>=2));
        addParameter(ip,'Kmax',8, @(z)isnumeric(z) && isscalar(z) && z>=2);
        addParameter(ip,'Distance','euclidean', @ischar);
        addParameter(ip,'Linkage','ward', @ischar);
        addParameter(ip,'Epsilon',[], @(z)isempty(z) || (isscalar(z) && z>0));
        addParameter(ip,'MinPts',[], @(z)isempty(z) || (isscalar(z) && z>=3));
        addParameter(ip,'UseUspace',false, @islogical);
        addParameter(ip,'CopulaFamilies', {'Gaussian','t','Clayton','Frank','Gumbel'}, @(c)iscellstr(c) || all(cellfun(@ischar,c)));
        addParameter(ip,'FitMixture2D',false, @islogical);
        addParameter(ip,'MixMaxComp',3, @(z)isscalar(z) && z>=2 && z<=8);
        addParameter(ip,'MixMaxIter',200, @(z)isscalar(z) && z>=20);
        addParameter(ip,'MixTol',1e-6, @(z)isnumeric(z) && z>0);
        addParameter(ip,'Verbose',true, @islogical);
        addParameter(ip,'Seed',[], @(z)isempty(z) || (isscalar(z) && isnumeric(z)));
        parse(ip, X, varargin{:});
        P = ip.Results;

        if ~isempty(P.Seed)
            try, rng(P.Seed,'twister'); end
        end

        n0 = size(X,1);
        X = X(~any(~isfinite(X),2),:);
        if size(X,1) < n0
            warn_('Dropped rows with NaN/Inf.');
        end
        U = makeUfromMarginals(X, P.Marginals);
        Z = X; if P.UseUspace, Z = U; end

        alg = lower(P.Algorithm);
        switch alg
            case 'kmeans'
                [clust,Kchosen,score] = run_kmeans(Z, P);
            case 'gmm'
                [clust,Kchosen,score] = run_gmm(Z, P);
            case 'hier'
                [clust,Kchosen,score] = run_hier(Z, P);
            case 'dbscan'
                [clust,Kchosen,score] = run_dbscan(Z, P);
            otherwise, error('Unknown algorithm.');
        end
        [clust, Kchosen, noiseIdx] = normalize_labels(clust, alg);
        if P.Verbose
            fprintf('[Cluster] Algorithm=%s | K=%d | n=%d | d=%d\n', alg, Kchosen, size(X,1), size(X,2));
            if isfield(score,'sil'), fprintf('  Silhouette=%.3f\n',score.sil); end
            if isfield(score,'gap'), fprintf('  Gap=%.3f\n',score.gap); end
            if isfield(score,'bic'), fprintf('  BIC=%.1f\n',score.bic); end
        end

        famAllow = lower(unique(P.CopulaFamilies));
        clusters = repmat(struct(), Kchosen, 1);
        for k=1:Kchosen
            idx = (clust==k);
            Xk  = X(idx,:); Uk = U(idx,:);
            clusters(k).idx = find(idx);
            clusters(k).size= sum(idx);
            clusters(k).X   = Xk;
            clusters(k).U   = Uk;
            fams = famAllow;
            if size(X,2)>2
                fams = intersect(fams, {'gaussian','t'});
                if isempty(fams), fams = {'gaussian','t'}; end
            end
            [fitK, critK] = fit_best_copula(Uk, fams);
            clusters(k).fit     = fitK;
            clusters(k).criteria= critK;
        end

        mixture2D = [];
        if P.FitMixture2D && size(U,2)==2
            mixture2D = em_mixture_copulas(U, famAllow, P.MixMaxComp, P.MixMaxIter, P.MixTol, P.Verbose);
        end

        plots = struct();
        try, plots.clusterScatter = plot_cluster_scatter(X, clust, noiseIdx); catch ME, warn_(ME.message); end
        try, plots.clusterUContours = plot_cluster_U_contours(U, clust, clusters); catch ME, warn_(ME.message); end

        CC = struct();
        CC.X = X; CC.U = U; CC.marginalName = P.Marginals;
        CC.alg = alg; CC.K = Kchosen; CC.clust = clust; CC.scores = score;
        CC.clusters = clusters; CC.noiseIdx = noiseIdx; CC.mixture2D = mixture2D; CC.plots = plots;

        % --- local helpers (inline module) ---
        function [labels,Kbest,score] = run_kmeans(Z, P)
            Kcand = decide_K_candidates(P.K, P.Kmax, Z);
            best = struct('K',[], 'sil',-Inf, 'gap',-Inf, 'labels',[]);
            for K = Kcand
                try
                    opts = statset('MaxIter',1000,'UseParallel',false,'Display','off');
                    [lab,~] = kmeans(Z, K, 'Replicates',10, 'Distance','sqeuclidean', 'MaxIter',1000, 'Options',opts, 'OnlinePhase','on');
                    sil = mean(silhouette_quick(Z, lab));
                    gap = gap_statistic(Z, lab);
                    scoreK = sil + 0.15*gap; %#ok<NASGU>
                    if sil > best.sil
                        best.K = K; best.sil = sil; best.gap = gap; best.labels = lab;
                    end
                catch
                end
            end
            labels = best.labels;
            if isempty(labels)
                labels = kmeans(Z, 2, 'Replicates',5);
                best.K = 2; best.sil = mean(silhouette_quick(Z, labels)); best.gap = gap_statistic(Z, labels);
            end
            Kbest = best.K;
            score = struct('sil',best.sil,'gap',best.gap);
        end
        function [labels,Kbest,score] = run_gmm(Z, P)
            Kcand = decide_K_candidates(P.K, P.Kmax, Z);
            best = struct('K',[], 'bic',Inf, 'labels',[]);
            for K = Kcand
                try
                    GM = fitgmdist(Z, K, 'CovarianceType','full','RegularizationValue',1e-6,'Replicates',10, 'Options',statset('MaxIter',1000));
                    [~,~,post] = cluster(GM, Z); [~,lab] = max(post,[],2);
                    bic = GM.BIC;
                    if bic < best.bic
                        best.K = K; best.bic = bic; best.labels = lab;
                    end
                catch
                end
            end
            labels = best.labels;
            if isempty(labels)
                GM = fitgmdist(Z, 2, 'CovarianceType','full','RegularizationValue',1e-6,'Replicates',5);
                [~,~,post] = cluster(GM, Z); [~,labels] = max(post,[],2);
                best.K = 2; best.bic = GM.BIC;
            end
            Kbest = best.K;
            score = struct('bic',best.bic);
        end
        function [labels,Kbest,score] = run_hier(Z, P)
            D = pdist(Z,'euclidean');
            L = linkage(D, P.Linkage);
            Kcand = decide_K_candidates(P.K, P.Kmax, Z);
            best = struct('K',[], 'sil',-Inf, 'labels',[]);
            for K = Kcand
                lab = cluster(L,'Maxclust',K);
                sil = mean(silhouette_quick(Z, lab));
                if sil > best.sil
                    best.K = K; best.sil = sil; best.labels = lab;
                end
            end
            labels = best.labels;
            if isempty(labels)
                labels = cluster(L,'Maxclust',2);
                best.K = 2; best.sil = mean(silhouette_quick(Z, labels));
            end
            Kbest = best.K; score = struct('sil',best.sil);
        end
        function [labels,Kbest,score] = run_dbscan(Z, P)
            eps = P.Epsilon; minPts = P.MinPts;
            if isempty(minPts), minPts = max(5, ceil(size(Z,2)*2)); end
            if isempty(eps), eps = estimate_eps(Z, minPts); end
            labels = dbscan(Z, eps, minPts, 'Distance','euclidean');
            Kbest = max(labels);
            score = struct('sil',NaN,'gap',NaN);
        end
        function Kcand = decide_K_candidates(K, Kmax, Z)
            if ~isempty(K)
                Kcand = K(:)';
            else
                Kcand = 2:min(Kmax, max(8, ceil(sqrt(size(Z,1))/2)));
                Kcand = unique(Kcand);
            end
        end
        function eps = estimate_eps(Z, minPts)
            D = pdist2(Z, Z);
            D(1:size(D,1)+1:end)=Inf;
            dk = zeros(size(Z,1),1);
            for ii=1:size(Z,1)
                di = sort(D(ii,:),'ascend');
                dk(ii) = di(min(minPts, numel(di)));
            end
            eps = prctile(dk, 80);
            eps = max(eps, 1e-6);
        end
        function silv = silhouette_quick(Z, lab)
            n = size(Z,1);
            if n>5000
                idx = randperm(n, 5000);
                silv = silhouette(Z(idx,:), lab(idx), 'Euclidean');
            else
                silv = silhouette(Z, lab, 'Euclidean');
            end
            silv = silv(:); silv = silv(isfinite(silv));
            if isempty(silv), silv = 0; end
        end
        function g = gap_statistic(Z, lab)
            K = max(lab); if K<2, g=0; return; end
            W = within_dispersion(Z, lab);
            bounds = [min(Z,[],1); max(Z,[],1)];
            Wref = zeros(10,1);
            for b=1:10
                Zr = rand_uniform_in_box(size(Z,1), bounds);
                labr = kmeans(Zr, K, 'Replicates',3, 'MaxIter',300, 'Display','off');
                Wref(b) = within_dispersion(Zr, labr);
            end
            g = mean(log(Wref)) - log(W);
        end
        function Zr = rand_uniform_in_box(n, bounds)
            lo = bounds(1,:); hi = bounds(2,:);
            Zr = rand(n, numel(lo)).*(hi-lo) + lo;
        end
        function W = within_dispersion(Z, lab)
            K = max(lab); W=0;
            for kk=1:K
                z = Z(lab==kk,:);
                if size(z,1)>1
                    c = mean(z,1);
                    W = W + sum(sum((z-c).^2));
                end
            end
        end
        function [lab, K, noiseIdx] = normalize_labels(lab, alg)
            noiseIdx = find(lab==0);
            if strcmpi(alg,'dbscan')
                K = max(lab);
            else
                K = max(lab);
            end
            nz = find(lab>0);
            if ~isempty(nz)
                [~,~,labn] = unique(lab(nz));
                lab2 = zeros(size(lab));
                lab2(nz) = labn;
                lab = lab2;
                K = max(lab);
            end
        end
        function [fitBest, critBest] = fit_best_copula(Uk, fams)
            n = size(Uk,1); d = size(Uk,2);
            best = struct('fam','','logL',-Inf,'AIC',Inf,'BIC',Inf,'fit',[]);
            for f = 1:numel(fams)
                fam = lower(fams{f});
                try
                    [fit, logL, kparams] = single_fit(Uk, fam);
                    AIC = -2*logL + 2*kparams;
                    BIC = -2*logL + kparams*log(n);
                    if BIC < best.BIC
                        best.fam = fam; best.logL = logL; best.AIC = AIC; best.BIC = BIC; best.fit = fit;
                    end
                catch
                end
            end
            fitBest = best.fit;
            critBest = rmfield(best, 'fit');
            if ~isempty(fitBest)
                fitBest.family = best.fam;
            end
            function [fit, logL, kparams] = single_fit(U, fam)
                switch fam
                    case 'gaussian'
                        R = copulafit('Gaussian', U);
                        R = nearestPD(R);
                        c = copulapdf('Gaussian', U, R);
                        logL = sum(log(max(c, realmin)));
                        kparams = nnz(triu(true(size(R)),1));
                        fit = struct('family','gaussian','R',R,'nu',[],'theta',[]);
                    case 't'
                        [R, nu] = copulafit('t', U);
                        R = nearestPD(R); nu = max(2,round(nu));
                        c = copulapdf('t', U, R, nu);
                        logL = sum(log(max(c, realmin)));
                        kparams = nnz(triu(true(size(R)),1)) + 1;
                        fit = struct('family','t','R',R,'nu',nu,'theta',[]);
                    case {'clayton','frank','gumbel'}
                        if size(U,2) ~= 2, error('Archimedean supports only 2D here.'); end
                        th = copulafit(ucfirst(fam), U);
                        c = copulapdf(ucfirst(fam), U, th);
                        logL = sum(log(max(c, realmin)));
                        kparams = 1;
                        fit = struct('family',fam,'R',[],'nu',[],'theta',th);
                    otherwise
                        error('Unknown family %s', fam);
                end
            end
        end
        function outmix = em_mixture_copulas(U, famAllow, Kmax, maxIter, tol, verbose)
            if size(U,2) ~= 2, outmix = []; return; end
            fams = lower(unique(famAllow)); fams = fams(:)';
            n = size(U,1);
            bestBIC = Inf; best = [];
            for K = 2:Kmax
                [lab,~] = kmeans(U, K, 'Replicates',5,'MaxIter',500,'Display','off');
                Rsp = zeros(n,K);
                for kk=1:K, Rsp(:,kk) = (lab==kk); end
                Rsp = bsxfun(@rdivide, Rsp, sum(Rsp,2)+eps);
                compFam = repmat({'gaussian'}, 1, K);
                for kk=2:2:K, compFam{kk} = 't'; end
                params = cell(1,K);
                ll_old = -Inf;
                for it=1:maxIter
                    for kk=1:K
                        wk = Rsp(:,kk);
                        params{kk} = weighted_fit(U, wk, fams);
                    end
                    pi_k = mean(Rsp,1);
                    C = zeros(n,K);
                    for kk=1:K
                        C(:,kk) = max(component_pdf(U, params{kk}), realmin).*max(pi_k(kk),1e-9);
                    end
                    sumC = sum(C,2); ll = sum(log(sumC));
                    Rsp = bsxfun(@rdivide, C, sumC + realmin);
                    if verbose && mod(it,25)==0, fprintf('  EM[K=%d] iter=%d ll=%.3f\n',K,it,ll); end
                    if abs(ll-ll_old) < tol*max(1,abs(ll_old)), break; end
                    ll_old = ll;
                end
                kparams = 0;
                for kk=1:K, kparams = kparams + num_params(params{kk}); end
                kparams = kparams + (K-1);
                BIC = -2*ll + kparams*log(n);
                if BIC < bestBIC
                    bestBIC = BIC;
                    best = struct('K',K,'params',{params},'pi',pi_k,'logL',ll,'BIC',BIC,'Rsp',Rsp);
                end
            end
            outmix = best;
            function p = component_pdf(U, prm)
                switch prm.family
                    case 'gaussian', p = copulapdf('Gaussian', U, prm.R);
                    case 't',        p = copulapdf('t', U, prm.R, prm.nu);
                    case 'clayton',  p = copulapdf('Clayton', U, prm.theta);
                    case 'frank',    p = copulapdf('Frank', U, prm.theta);
                    case 'gumbel',   p = copulapdf('Gumbel', U, prm.theta);
                    otherwise,       p = ones(size(U,1),1);
                end
                p(~isfinite(p)) = realmin;
            end
            function prm = weighted_fit(U, w, fams)
                w = w(:); w = w./(sum(w)+eps);
                famsTry = fams;
                bestw = struct('family','','logL',-Inf,'prm',[]);
                for ff=1:numel(famsTry)
                    fam = famsTry{ff};
                    try
                        switch fam
                            case 'gaussian'
                                R = copulafit('Gaussian', U, 'Weights', w);
                                c = copulapdf('Gaussian', U, R); ll = sum(w.*log(max(c,realmin)));
                                prm0 = struct('family','gaussian','R',nearestPD(R),'nu',[],'theta',[]);
                                kparams = 1; %#ok<NASGU>
                            case 't'
                                [R,nu] = copulafit('t', U, 'Weights', w);
                                R = nearestPD(R); nu = max(2,round(nu));
                                c = copulapdf('t', U, R, nu); ll = sum(w.*log(max(c,realmin)));
                                prm0 = struct('family','t','R',R,'nu',nu,'theta',[]);
                            case {'clayton','frank','gumbel'}
                                th = copulafit(ucfirst(fam), U, 'Weights', w);
                                c = copulapdf(ucfirst(fam), U, th); ll = sum(w.*log(max(c,realmin)));
                                prm0 = struct('family',fam,'R',[],'nu',[],'theta',th);
                            otherwise
                                continue
                        end
                        if ll > bestw.logL
                            bestw.family = fam; bestw.logL = ll; bestw.prm = prm0;
                        end
                    catch
                    end
                end
                prm = bestw.prm;
                if isempty(prm), prm = struct('family','gaussian','R',eye(2),'nu',[],'theta',[]); end
            end
            function k = num_params(prm)
                switch prm.family
                    case 'gaussian', k = 1;
                    case 't',        k = 2;
                    otherwise,       k = 1;
                end
            end
        end
        function h = plot_cluster_scatter(X, clust, noiseIdx)
            h = figure('Name','Clusters (raw space)','Color','w');
            d = size(X,2);
            if d==2
                gscatter(X(:,1),X(:,2),clust); grid on;
                xlabel('X1'); ylabel('X2'); title('Clusters (raw)');
                if ~isempty(noiseIdx), hold on; plot(X(noiseIdx,1),X(noiseIdx,2),'kx'); legend('show'); end
            else
                [coeff,score,~,~,expl] = pca(X,'NumComponents',3); %#ok<ASGLU>
                scatter3(score(:,1),score(:,2),score(:,3),24,clust,'filled'); grid on;
                xlabel(sprintf('PC1 (%.1f%%)',expl(1))); ylabel(sprintf('PC2 (%.1f%%)',expl(2))); zlabel(sprintf('PC3 (%.1f%%)',expl(3)));
                title('Clusters on PCA space');
                colorbar;
            end
        end
        function h = plot_cluster_U_contours(U, clust, clusters)
            if size(U,2) ~= 2, h = []; return; end
            h = figure('Name','Cluster-wise copula contours (U-space)','Color','w');
            K = max(clust);
            tlay = tiledlayout(ceil(K/2),2, 'TileSpacing','compact','Padding','compact');
            for kk=1:K
                nexttile;
                Uk = U(clust==kk,:);
                if isempty(Uk), text(0.5,0.5,sprintf('Cluster %d empty',kk)); axis off; continue; end
                g = linspace(0.01,0.99,80); [G1,G2] = meshgrid(g,g);
                fit = clusters(kk).fit;
                switch lower(fit.family)
                    case 'gaussian', Cpdf = copulapdf('Gaussian',[G1(:) G2(:)], fit.R);
                    case 't',        Cpdf = copulapdf('t',[G1(:) G2(:)], fit.R, fit.nu);
                    case 'clayton',  Cpdf = copulapdf('Clayton',[G1(:) G2(:)], fit.theta);
                    case 'frank',    Cpdf = copulapdf('Frank',[G1(:) G2(:)], fit.theta);
                    case 'gumbel',   Cpdf = copulapdf('Gumbel',[G1(:) G2(:)], fit.theta);
                    otherwise,       Cpdf = ones(numel(G1),1);
                end
                Cpdf = reshape(Cpdf, size(G1));
                contourf(g,g,Cpdf,12,'LineStyle','none'); colorbar;
                hold on; scatter(Uk(:,1),Uk(:,2),8,'k','filled','MarkerFaceAlpha',0.2);
                xlim([0 1]); ylim([0 1]); title(sprintf('Cluster %d (%s)',kk, ucfirst(fit.family)));
            end
            title(tlay,'Copula PDF per cluster (U-space)');
        end
        function s = ucfirst(fm), s = lower(fm); if ~isempty(s), s(1)=upper(s(1)); end, end
        function warn_(msg)
            try, warning('CC:warn','%s',msg); catch, fprintf(2,'[WARN] %s\n',msg); end
        end
    end

%% ---- Cluster Validity Dialog ----
function clusterValidityDialog(src)
    try
        CC = []; if isfield(userData,'Cluster') && isfield(userData.Cluster,'results'), CC = userData.Cluster.results; end
        if isempty(CC) || ~isfield(CC,'labels') || isempty(CC.labels)
            showAlert(src,'Run "Cluster & Fit Copulas" first.','Cluster Validity'); return;
        end
        cols = getSelectedColumns(); X = userData.Table{:,cols}; marg = ddlMarg.Value;
        useU = true; if isfield(CC,'useU'), useU = CC.useU; end
        if useU, Z = makeUfromMarginals(X, marg); else, Z = X; end
        labels = CC.labels(:);
        stats = compute_cluster_validity(Z, labels);
        fig = uifigure('Name','Cluster Validity','Position',[340 340 540 280]);
        uitable(fig,'Data', {stats.Silhouette, stats.Gap, stats.CalinskiHarabasz, stats.DaviesBouldin},...
            'ColumnName',{'Silhouette','Gap','Calinski–Harabasz','Davies–Bouldin'},'Position',[20 160 500 50]);
        % Bar plot
        figure('Name','Validity Indices'); vals = [stats.Silhouette, stats.Gap, stats.CalinskiHarabasz, 1/max(stats.DaviesBouldin,eps)];
        bar(vals); set(gca,'XTickLabel',{'Sil','Gap','CH','1/DB'}); grid on; title('Cluster Validity');
    catch ME, showAlert(src,ME.message,'Validity Error'); end
end

function stats = compute_cluster_validity(Z, labels)
    labels = labels(:); K = numel(unique(labels));
    N = size(Z,1);
    sil = NaN; gapv = NaN; ch = NaN; db = NaN;
    try
        svals = silhouette(Z, labels);
        sil = mean(svals,'omitnan');
    catch, sil = NaN; end
    try
        % Simple gap approximation via log(W_k)
        Wk = 0;
        for k=1:K
            idx = labels==k; Zk = Z(idx,:);
            if size(Zk,1)>=2
                Ck = pdist2(Zk, mean(Zk,1)); Wk = Wk + sum(Ck.^2);
            end
        end
        gapv = -log(max(Wk,eps));
    catch, gapv = NaN; end
    try
        ch = evalclusters(Z,labels,'CalinskiHarabasz'); ch = ch.CriterionValues; ch = ch(end);
        if iscell(ch), ch = ch{end}; end
        ch = double(ch);
    catch, ch = NaN; end
    try
        db = evalclusters(Z,labels,'DaviesBouldin'); db = db.CriterionValues; db = db(end);
        if iscell(db), db = db{end}; end
        db = double(db);
    catch, db = NaN; end
    stats = struct('Silhouette',sil,'Gap',gapv,'CalinskiHarabasz',ch,'DaviesBouldin',db);
end

%% ---- Silhouette per cluster ----
function plotSilhouettePerCluster(src)
    try
        CC = []; if isfield(userData,'Cluster') && isfield(userData.Cluster,'results'), CC = userData.Cluster.results; end
        if isempty(CC) || ~isfield(CC,'labels') || isempty(CC.labels)
            showAlert(src,'Run "Cluster & Fit Copulas" first.','Silhouette'); return;
        end
        cols = getSelectedColumns(); X = userData.Table{:,cols}; marg = ddlMarg.Value;
        Z = X; if CC.useU, Z = makeUfromMarginals(X, marg); end
        [~,axh] = silhouette(Z, CC.labels);
        title(axh, 'Silhouette by point'); grid(axh,'on');
    catch ME, showAlert(src,ME.message,'Silhouette Error'); end
end

%% ---- Consensus Clustering ----
function consensusClusteringDialog(src)
    try
        CC = []; if isfield(userData,'Cluster') && isfield(userData.Cluster,'results'), CC = userData.Cluster.results; end
        if isempty(CC) || ~isfield(CC,'labels') || isempty(CC.labels)
            showAlert(src,'Run "Cluster & Fit Copulas" first.','Consensus'); return;
        end
        B = 100; % fixed for now
        cols = getSelectedColumns(); X = userData.Table{:,cols}; marg = ddlMarg.Value;
        Z = X; if CC.useU, Z = makeUfromMarginals(X, marg); end
        K = CC.K; alg = CC.algorithm;
        [S,labels_cons] = consensus_clustering(Z, K, B, alg);
        userData.Cluster.consensus = struct('S',S,'labels',labels_cons,'B',B);
        figure('Name','Consensus Co-association'); imagesc(S); colorbar; title('Consensus S (co-association)');
        showAlert(src, 'Consensus clustering computed (B=100).','Consensus');
    catch ME, showAlert(src,ME.message,'Consensus Error'); end
end

function [S, labels_cons] = consensus_clustering(Z, K, B, alg)
    n = size(Z,1); S = zeros(n);
    for b=1:B
        idx = randsample(n, n, true);
        Zb = Z(idx,:);
        lb = run_single_clustering(Zb, K, alg);
        % map back votes
        for i=1:n
            for j=i:n
                same = (lb(i)==lb(j));
                S(i,j) = S(i,j) + same;
                if i~=j, S(j,i)=S(i,j); end
            end
        end
    end
    S = S / B;
    % cluster S into K clusters
    D = 1 - S;
    labels_cons = spectralClusterWrapper(D, K);
end

function labels = run_single_clustering(Z, K, alg)
    switch lower(alg)
        case 'kmeans'
            labels = kmeans(Z, K, 'Replicates',3,'MaxIter',300,'OnlinePhase','off');
        case 'gmm'
            GM = fitgmdist(Z, K, 'RegularizationValue',1e-6,'Options',statset('MaxIter',500));
            labels = cluster(GM, Z);
        case 'hier'
            L = linkage(Z,'ward'); labels = cluster(L,'maxclust',K);
        case 'dbscan'
            % Fallback: approximate using kmeans when eps not provided
            labels = kmeans(Z, K, 'Replicates',3,'MaxIter',300,'OnlinePhase','off');
        otherwise
            labels = kmeans(Z, K, 'Replicates',3);
    end
end

%% ---- Spectral clustering with copula-based distance ----
function spectralClusteringDialog(src)
    try
        cols = getSelectedColumns(); if numel(cols)<2, showAlert(src,'Select ≥2 columns.','Spectral'); return; end
        X = userData.Table{:,cols}; marg = ddlMarg.Value; U = makeUfromMarginals(X, marg);
        % Distance: empirical copula discrepancy (fast approx via ranks)
        R = tiedrank(U) / (size(U,1)+1);
        % Build affinity via RBF on rank space
        D = squareform(pdist(R,'euclidean'));
        sigma = median(D(D>0)); if ~isfinite(sigma) || sigma<=0, sigma = 1; end
        A = exp(-(D.^2)/(2*sigma^2));
        L = diag(sum(A,2)) - A;
        K = max(2, min(10, numel(unique(userData.Cluster.results.labels))));
        labels = spectralClusterWrapper(L, K);
        figure('Name','Spectral clustering (copula distance)'); gscatter(U(:,1),U(:,2),labels); grid on; title('Spectral on U (RBF ranks)');
    catch ME, showAlert(src,ME.message,'Spectral Error'); end
end

function labels = spectralClusterWrapper(DorL, K)
    % Accepts distance matrix D or Laplacian L; detect via symmetry zero-diag
    L = DorL;
    if all(diag(DorL)==0) && all(all(DorL==DorL.')) && min(DorL(DorL>0))>=0 %#ok<AND2>
        % looks like a distance matrix -> build RBF affinity
        D = DorL; sigma = median(D(D>0)); if ~isfinite(sigma) || sigma<=0, sigma=1; end
        A = exp(-(D.^2)/(2*sigma^2));
        L = diag(sum(A,2)) - A;
    end
    [V,~] = eigs(L, K, 'SM'); Y = bsxfun(@rdivide, V, sqrt(sum(V.^2,2))+eps);
    labels = kmeans(Y, K, 'Replicates',10,'MaxIter',300,'OnlinePhase','off');
end

%% ---- Kernel k-means (HSIC) ----
function kernelKMeansDialog(src)
    try
        CC = []; if isfield(userData,'Cluster') && isfield(userData.Cluster,'results'), CC = userData.Cluster.results; end
        cols = getSelectedColumns(); X = userData.Table{:,cols}; marg = ddlMarg.Value; U = makeUfromMarginals(X, marg);
        if isempty(CC) || ~isfield(CC,'K'), K = max(2, round(sqrt(size(U,1))/5)); else, K = CC.K; end
        sigma = median(pdist(U)); if ~isfinite(sigma) || sigma<=0, sigma=1; end
        labels = kernel_kmeans(U, K, sigma);
        figure('Name','Kernel k-means (HSIC)'); gscatter(U(:,1),U(:,2),labels); grid on; title(sprintf('Kernel k-means (K=%d, sigma=%.3f)',K,sigma));
    catch ME, showAlert(src,ME.message,'Kernel k-means Error'); end
end

function labels = kernel_kmeans(U, K, sigma)
    n = size(U,1);
    D = squareform(pdist(U,'euclidean'));
    Kmat = exp(-(D.^2)/(2*sigma^2));
    labels = randi(K,n,1);
    for it=1:50
        H = zeros(n,K);
        for k=1:K, H(:,k) = (labels==k); end
        Nk = sum(H,1) + eps;
        dist = zeros(n,K);
        for k=1:K
            hk = H(:,k)/Nk(k);
            dist(:,k) = diag(Kmat) - 2*(Kmat*hk) + (hk.'*Kmat*hk);
        end
        [~,labels] = min(dist,[],2);
    end
end

%% ---- Temporal regimes (rolling clustering) ----
function temporalRegimesDialog(src)
    try
        cols = getSelectedColumns(); if numel(cols)<2, showAlert(src,'Select ≥2 columns.','Temporal'); return; end
        X = userData.Table{:,cols}; marg = ddlMarg.Value; U = makeUfromMarginals(X, marg);
        n = size(U,1);
        W = max(60, min(250, floor(n/5))); step = max(5, round(W/10));
        K = max(2, min(6, isfield(userData,'Cluster') * userData.Cluster.results.K + 0));
        labels_time = rolling_clustering(U, W, step, K);
        figure('Name','Temporal Regimes'); imagesc(labels_time.'); colorbar; xlabel('window idx'); ylabel('cluster'); title('Temporal cluster ids');
    catch ME, showAlert(src,ME.message,'Temporal Error'); end
end

function labels_time = rolling_clustering(U, W, step, K)
    n = size(U,1); win = 1:W; t = 1; labels_time = [];
    while win(end)<=n
        Uw = U(win,:);
        lb = kmeans(Uw, K, 'Replicates',3,'MaxIter',300,'OnlinePhase','off');
        labels_time(end+1,1) = mode(lb); %#ok<AGROW>
        win = win + step; t = t+1;
    end
end

%% ---- Tail Heatmap per cluster ----
function tailHeatmapPerCluster(src)
    try
        CC = []; if isfield(userData,'Cluster') && isfield(userData.Cluster,'results'), CC = userData.Cluster.results; end
        if isempty(CC) || ~isfield(CC,'clusters') || isempty(CC.clusters)
            showAlert(src,'Run "Cluster & Fit Copulas" first.','Tail Heatmap'); return;
        end
        cols = getSelectedColumns(); d = numel(cols);
        for k = 1:CC.K
            fitk = CC.clusters(k).fit; if ~isfield(fitk,'family'), continue; end
            if any(strcmpi(fitk.family,{'Gaussian','t'}))
                figure('Name',sprintf('Tail heatmap - cluster %d',k));
                R = fitk.R; imagesc(R); colorbar; title(sprintf('R (proxy tails) - %s', fitk.family));
            end
        end
    catch ME, showAlert(src,ME.message,'Tail Heatmap Error'); end
end

%% ---- Cluster Profile ----
function clusterProfileDialog(src)
    try
        CC = []; if isfield(userData,'Cluster') && isfield(userData.Cluster,'results'), CC = userData.Cluster.results; end
        if isempty(CC) || ~isfield(CC,'labels') || isempty(CC.labels)
            showAlert(src,'Run "Cluster & Fit Copulas" first.','Cluster Profile'); return;
        end
        cols = getSelectedColumns(); X = userData.Table{:,cols};
        K = CC.K; out = cell(K, 4);
        for k=1:K
            idx = CC.labels==k; Xk = X(idx,:);
            med = median(Xk,1,'omitnan'); iqrk = iqr(Xk);
            fam = CC.clusters(k).fit.family;
            bic = CC.clusters(k).criteria.BIC;
            out{k,1} = k; out{k,2} = sum(idx); out{k,3} = fam; out{k,4} = bic;
        end
        fig = uifigure('Name','Cluster Profile','Position',[320 320 520 240]);
        uitable(fig,'Data',out,'ColumnName',{'Cluster','n','Family','BIC'},'Position',[20 20 480 200]);
    catch ME, showAlert(src,ME.message,'Profile Error'); end
end

%% ---- Auto K (grid search over K) ----
function autoKGridDialog(src)
    try
        cols = getSelectedColumns(); if numel(cols)<2, showAlert(src,'Select ≥2 columns.','Auto K'); return; end
        X = userData.Table{:,cols}; marg = ddlMarg.Value;
        useU = true; if isfield(userData,'Cluster') && isfield(userData.Cluster,'results') && isfield(userData.Cluster.results,'useU'), useU = userData.Cluster.results.useU; end
        Z = X; if useU, Z = makeUfromMarginals(X, marg); end
        Kmax = min(12, max(3, round(sqrt(size(Z,1))/2)));
        Kgrid = 2:Kmax;
        silv = NaN(size(Kgrid)); chv = NaN(size(Kgrid)); dbv = NaN(size(Kgrid));
        for i=1:numel(Kgrid)
            K = Kgrid(i);
            lb = kmeans(Z, K, 'Replicates',5,'MaxIter',300,'OnlinePhase','off');
            try, silv(i) = mean(silhouette(Z, lb),'omitnan'); end
            try, ev = evalclusters(Z,lb,'CalinskiHarabasz'); chv(i) = ev.CriterionValues(end); end
            try, ev2 = evalclusters(Z,lb,'DaviesBouldin'); dbv(i) = ev2.CriterionValues(end); end
        end
        figure('Name','Auto K diagnostics'); 
        plot(Kgrid, silv,'-o'); hold on; plot(Kgrid, normalize_pos(chv),'-s'); plot(Kgrid, 1./max(dbv,eps),'-d'); grid on;
        xlabel('K'); ylabel('score (scaled)'); legend('Silhouette','CH (scaled)','1/DB'); title('Auto K grid');
        [~,iSil] = max(silv); [~,iCH] = max(chv); [~,iDB] = max(1./max(dbv,eps));
        msg = sprintf('Auto-K suggestions:\n  Silhouette → K=%d\n  CH → K=%d\n  1/DB → K=%d', Kgrid(iSil), Kgrid(iCH), Kgrid(iDB));
        showAlert(src, msg, 'Auto K');
    catch ME, showAlert(src,ME.message,'Auto K Error'); end
end

function y = normalize_pos(x)
    x = x(:); if all(~isfinite(x)), y = x; return; end
    m = min(x(isfinite(x))); M = max(x(isfinite(x))); if M<=m, y = ones(size(x)); else, y = (x-m)/(M-m); end
end

%% ---- DBSCAN epsilon elbow (k-dist curve) ----
function dbscanElbowDialog(src)
    try
        cols = getSelectedColumns(); if numel(cols)<2, showAlert(src,'Select ≥2 columns.','DBSCAN ε'); return; end
        X = userData.Table{:,cols}; marg = ddlMarg.Value; Z = makeUfromMarginals(X, marg);
        MinPts = max(5, round(log(size(Z,1)))); % heuristic
        D = squareform(pdist(Z,'euclidean'));
        Dsort = sort(D,2);
        Dk = Dsort(:,MinPts+1);
        [~,idx] = sort(Dk,'ascend'); y = Dk(idx); x = (1:numel(y))';
        figure('Name','DBSCAN k-distance curve'); plot(x,y,'-'); grid on; title(sprintf('k-distance (k=%d)',MinPts));
        kbest = elbow_maxcurv(y);
        yel = y(kbest); showAlert(src, sprintf('\\epsilon ≈ %.3f (knee at %d-th point)', yel, kbest),'DBSCAN ε elbow');
        hold on; plot(kbest, yel, 'o'); 
    catch ME, showAlert(src,ME.message,'DBSCAN ε Error'); end
end

function kidx = elbow_maxcurv(y)
    y = y(:); n = numel(y); if n<5, kidx = round(n/2); return; end
    x = (1:n)'; y1 = [y(2:end); y(end)]; y0 = [y(1); y(1:end-1)];
    curv = abs(y1 - 2*y + y0);
    [~,kidx] = max(curv(2:end-1)); kidx = kidx + 1;
end

%% ---- Outlier prune in U-space ----
function outlierPruneUDialog(src)
    try
        cols = getSelectedColumns(); X = userData.Table{:,cols}; marg = ddlMarg.Value; U = makeUfromMarginals(X, marg);
        k = max(5, round(log(size(U,1))));
        D = squareform(pdist(U,'euclidean')); D(D==0)=eps;
        Dsort = sort(D,2);
        lrd = k ./ sum(Dsort(:,2:k+1),2);
        lof = (D * (1./lrd)) ./ (k * (1./lrd));
        thr = quantile(lof, 0.98);
        idxOut = find(lof>thr);
        userData.Cluster.outliers = idxOut;
        figure('Name','Outlier scan (U)'); scatter(U(:,1),U(:,2),12,[.6 .6 .6],'filled'); hold on;
        if ~isempty(idxOut), scatter(U(idxOut,1),U(idxOut,2),20,'r','filled'); end
        grid on; title(sprintf('Outliers in U (k=%d, thr=98%%): %d found',k,numel(idxOut)));
    catch ME, showAlert(src,ME.message,'Outlier Error'); end
end

%% ---- Per-cluster GOF (multi-tests on Rosenblatt W) ----

function perClusterGOFDialog(src)
    % Per-cluster GOF με σωστό Rosenblatt (2D) και ανθεκτικό fallback labels

    try
        % === 0) Πάρε αποτελέσματα clustering
        CC = [];
        if isfield(userData,'Cluster') && isfield(userData.Cluster,'results')
            CC = userData.Cluster.results;
        end
        if isempty(CC)
            showAlert(src,'Run "Cluster & Fit Copulas" first.','Per-cluster GOF'); 
            return;
        end

        % === 1) Πάρε labels με fallback (labels ή clust)
        lab = [];
        if isfield(CC,'labels') && ~isempty(CC.labels)
            lab = CC.labels(:);
        elseif isfield(CC,'clust') && ~isempty(CC.clust)
            lab = CC.clust(:);
        end
        if isempty(lab)
            showAlert(src,'Run "Cluster & Fit Copulas" first.','Per-cluster GOF'); 
            return;
        end

        % === 2) Χρειάζονται ακριβώς 2 στήλες (2D)
        cols = getSelectedColumns();
        if numel(cols)~=2
            showAlert(src,'Select exactly 2 columns (2D).','Per-cluster GOF');
            return;
        end

        % Τα δεδομένα + U-space με τα τρέχοντα marginals
        X    = userData.Table{:,cols};
        marg = ddlMarg.Value;
        U    = makeUfromMarginals(X, marg);

        % Αριθμός clusters
        K = max(lab);
        out = cell(K, 6);

        % Helper: archimedean proper name ('clayton' -> 'Clayton')
        archName = @(s) [upper(s(1)) lower(s(2:end))];

        % === 3) Loop ανά cluster
        for k = 1:K
            idx = (lab==k);
            Uk  = U(idx,:);
            nk  = sum(idx);

            % guard για μικρά clusters
            if nk < 30
                out(k,:) = {k, nk, NaN, NaN, NaN, NaN};
                continue;
            end

            % Πάρε τη fitted οικογένεια του cluster από το CC
            fam = 'gaussian';
            if isfield(CC,'clusters') && numel(CC.clusters) >= k ...
               && isfield(CC.clusters(k),'fit') && isfield(CC.clusters(k).fit,'family') ...
               && ~isempty(CC.clusters(k).fit.family)
                fam = lower(CC.clusters(k).fit.family);
            end

            % === 4) Rosenblatt transform με finite differences
            epsFD = 1e-4;
            U1 = Uk(:,1);
            U2 = Uk(:,2);

            % W(:,1)=U1, W(:,2)=∂C/∂u1(u1,u2)
            W = nan(size(Uk));
            W(:,1) = U1;

            switch fam
                case 'gaussian'
                    R = copulafit('Gaussian', Uk);
                    for i = 1:nk
                        u1 = U1(i); u2 = U2(i);
                        Cplus  = copulacdf('Gaussian', [min(u1+epsFD,1) u2], R);
                        Cminus = copulacdf('Gaussian', [max(u1-epsFD,0) u2], R);
                        V2 = (Cplus - Cminus) / (2*epsFD);
                        W(i,2) = max(0, min(1, V2));
                    end

                case 't'
                    [R,nu] = copulafit('t', Uk);
                    nu = max(2, round(nu));
                    for i = 1:nk
                        u1 = U1(i); u2 = U2(i);
                        Cplus  = copulacdf('t', [min(u1+epsFD,1) u2], R, nu);
                        Cminus = copulacdf('t', [max(u1-epsFD,0) u2], R, nu);
                        V2 = (Cplus - Cminus) / (2*epsFD);
                        W(i,2) = max(0, min(1, V2));
                    end

                case {'clayton','frank','gumbel'}
                    famMAT = archName(fam);                        % 'Clayton','Frank','Gumbel'
                    th = copulafit(famMAT, Uk);
                    for i = 1:nk
                        u1 = U1(i); u2 = U2(i);
                        Cplus  = copulacdf(famMAT, [min(u1+epsFD,1) u2], th);
                        Cminus = copulacdf(famMAT, [max(u1-epsFD,0) u2], th);
                        V2 = (Cplus - Cminus) / (2*epsFD);
                        W(i,2) = max(0, min(1, V2));
                    end

                otherwise
                    % άγνωστη ετικέτα -> γύρνα ως Gaussian fallback
                    R = copulafit('Gaussian', Uk);
                    for i = 1:nk
                        u1 = U1(i); u2 = U2(i);
                        Cplus  = copulacdf('Gaussian', [min(u1+epsFD,1) u2], R);
                        Cminus = copulacdf('Gaussian', [max(u1-epsFD,0) u2], R);
                        V2 = (Cplus - Cminus) / (2*epsFD);
                        W(i,2) = max(0, min(1, V2));
                    end
            end

            % === 5) GOF tests: U1 και V2 ~ U(0,1)
            [pKS,  pCvM,  pAD ] = gof_unif_tests(W(:,1));
            [pKS2, pCvM2, pAD2] = gof_unif_tests(W(:,2));
            pMin = min([pKS pCvM pAD pKS2 pCvM2 pAD2]);

            out(k,:) = {k, nk, pKS2, pCvM2, pAD2, pMin};
        end

        % === 6) Πίνακας αποτελεσμάτων
        fig = uifigure('Name','Per-cluster GOF','Position',[340 340 640 260]);
        uitable(fig, ...
            'Data', out, ...
            'ColumnName', {'Cluster','n','KS(V2)','CvM(V2)','AD(V2)','min p'}, ...
            'Position', [20 20 600 220]);

    catch ME
        showAlert(src, ME.message, 'Per-cluster GOF Error');
    end
end


function [pKS,pCvM,pAD] = gof_unif_tests(x)
    x = x(:); x = x(isfinite(x) & x>=0 & x<=1);
    if numel(x)<20, pKS=NaN; pCvM=NaN; pAD=NaN; return; end
    [~,pKS] = kstest(x, [sort(x) ( (1:numel(x))'/numel(x) )]);
    [pCvM, pAD] = cvm_ad_pvalues(x);
end

function [pCvM, pAD] = cvm_ad_pvalues(x)
    x = sort(x(:)); n = numel(x);
    ui = x; ti = (2*(1:n)'-1)./(2*n);
    W2 = (1/(12*n)) + sum( (ui - ti).^2 );
    A2 = -n - mean((2*(1:n)'-1).*log(ui) + (2*n+1-2*(1:n)').*log(1-ui));
    pCvM = exp(-4*W2); pAD = exp(-A2/2);
    pCvM = min(max(pCvM,0),1); pAD = min(max(pAD,0),1);
end

%% ---- Feature Importance (permutation on validity) ----
function featureImportanceDialog(src)
    try
        CC = []; 
if isfield(userData,'Cluster') && isfield(userData.Cluster,'results')
    CC = userData.Cluster.results;
end

% fallback για labels / algorithm / useU
lab = [];
if isfield(CC,'labels') && ~isempty(CC.labels)
    lab = CC.labels;
elseif isfield(CC,'clust') && ~isempty(CC.clust)
    lab = CC.clust;
end
if isempty(lab)
    showAlert(src,'Run "Cluster & Fit Copulas" first.','Importance'); 
    return;
end

cols = getSelectedColumns(); 
X = userData.Table{:,cols}; 
marg = ddlMarg.Value;

useU = false;
if isfield(CC,'useU'), useU = logical(CC.useU); end
Z = X; 
if useU
    Z = makeUfromMarginals(X, marg);
end

base_labels = lab(:);

        base_sil = mean(silhouette(Z, base_labels),'omitnan');
        imp = zeros(1,size(Z,2));
        for j=1:size(Z,2)
            Zp = Z; Zp(:,j) = Zp(randperm(size(Z,1)), j);
            sil_j = mean(silhouette(Zp, base_labels),'omitnan');
            imp(j) = max(0, base_sil - sil_j);
        end
        figure('Name','Feature Importance (permutation)'); bar(imp); grid on; 
        set(gca,'XTick',1:numel(cols),'XTickLabel',cols); title('Drop in Silhouette when permuted (higher = more important)');
    catch ME, showAlert(src,ME.message,'Importance Error'); end
end

%% ---- Export Cluster Report ----
function exportClusterReport(src)
    try
        CC = []; if isfield(userData,'Cluster') && isfield(userData.Cluster,'results'), CC = userData.Cluster.results; end
        if isempty(CC) || ~isfield(CC,'clusters') || isempty(CC.clusters)
            showAlert(src,'Run "Cluster & Fit Copulas" first.','Export Report'); return;
        end
        [file,path] = uiputfile('cluster_report.txt','Save Cluster Report'); if isequal(file,0), return; end
        fid = fopen(fullfile(path,file),'w');
        fprintf(fid,'Cluster Report\n'); fprintf(fid,'Alg=%s | K=%d | UseU=%d\n', CC.algorithm, CC.K, CC.useU);
        for k=1:CC.K
            cfit = CC.clusters(k).fit; bic = CC.clusters(k).criteria.BIC; n = sum(CC.labels==k);
            if isfield(cfit,'R') && ~isempty(cfit.R)
                fam = cfit.family; fprintf(fid,'Cluster %d: %s | n=%d | BIC=%.3f\n',k,fam,n,bic);
                fprintf(fid,'R=\n%s\n', evalc('disp(cfit.R)'));
                if isfield(cfit,'nu') && ~isempty(cfit.nu), fprintf(fid,'nu=%.3f\n', cfit.nu); end
            else
                fprintf(fid,'Cluster %d: %s | n=%d | BIC=%.3f\n',k,cfit.family,n,bic);
                if isfield(cfit,'theta') && ~isempty(cfit.theta), fprintf(fid,'theta=%.3f\n', cfit.theta); end
            end
            fprintf(fid,'---\n');
        end
        fclose(fid);
        showAlert(src,'Cluster report saved.','Export Cluster Report');
    catch ME, showAlert(src,ME.message,'Export Report Error'); end
end

%% ---- Spectral (Copula CvM distance) ----
function spectralCopulaCvMDialog(src)
    try
        cols = getSelectedColumns(); if numel(cols)<2, showAlert(src,'Select ≥2 columns.','Spectral CvM'); return; end
        X = userData.Table{:,cols}; marg = ddlMarg.Value; U = makeUfromMarginals(X, marg);
        n = size(U,1);
        D = squareform(pdist(tiedrank(U)/(n+1),'euclidean'));
        labels = spectralClusterWrapper(D, max(2, min(8, round(sqrt(n)/2))));
        figure('Name','Spectral (Copula CvM proxy)'); gscatter(U(:,1),U(:,2),labels); grid on; title('Spectral via copula distance proxy');
    catch ME, showAlert(src,ME.message,'Spectral CvM Error'); end
end


end
