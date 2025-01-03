function [T, parglmo] = parglm(X, F, varargin)

% MDSA EDIT: TAKING  TO ACCOUNT FOR COMPLEX NUMBERS
% Parallel General Linear Model to obtain multivariate factor and interaction 
% matrices in a crossed experimental design and permutation test for multivariate 
% statistical significance. Missing data is considered.
%
% Related routines: asca, apca, parglmVS, parglmMC, create_design
%
% T = parglm(X, F)   % minimum call
% [T, parglmo] = parglm(X, F,'Model',model,'Preprocessing',prep,'Permutations',n_perm,'Ts',ts,'Ordinal',ordinal,'Fmtc',fmtc,'Coding',coding,'Nested',nested)   % complete call
%
%
% INPUTS:
%
% X: [NxM] billinear data set for model fitting, where each row is a
% measurement, each column a variable
%
% F: [NxF] design matrix, cell or array, where columns correspond to 
% factors and rows to levels.
%
% Optional INPUTS (parameters):
%
% 'Model': This paremeter is similar to 'model' of anovan. It could be:
%       'linear': only main effects are provided (by default)
%       'interaction': two order interactions are provided
%       'full': all potential interactions are provided
%       [1x1]: maximum order of interactions considered
%       [ix2]: array with two order interactions
%       cell: with each element a vector of factors
%
% 'Preprocessing': [1x1] preprocesing:
%       0: no preprocessing 
%       1: mean-centering 
%       2: auto-scaling (default)  
%
% 'Permutations': [1x1] number of permutations (1000 by default)
%
% 'Ts': [1x1] Use SSQ (0) or the F-value (otherwise, by default) as test statistic  
%       0: Sum-of-squares of the factor/interaction
%       1: F-ratio of the SS of the factor/interaction divided by the SS of 
%       the residuals (by default)
%       2: F-ratio following the factors/interactions hierarchy (only for
%       random models or unconstrained mixed model, see Montogomery)
%
% 'Ordinal': [1xF] whether factors are nominal or ordinal
%       0: nominal (default)
%       1: ordinal
% 
% 'Fmtc': [1x1] correct for multiple-tesis when multifactorial (multi-way)
% analysis
%       0: do not correct (default)
%       1: Bonferroni 
%       2: Holm step-up or Hochberg step-down
%       3: Benjamini-Hochberg step-down (FDR)
%       4: Q-value from Benjamini-Hochberg step-down
%
% 'Coding': [1xF] type of coding of factors
%       0: sum/deviation coding (default)
%       1: reference coding (reference is the last level)
%
% 'Nested': [nx2] pairs of neted factors, e.g., if factor 2 is nested in 1,
%   and 3 in 2, then nested = [1 2; 2 3]
%
%
% OUTPUTS:
%
% T (table): ANOVA-like output table
%
% parglmo (structure): structure with the factor and interaction
% matrices, p-values (corrected, depending on fmtc) and explained variance 
%
%
% EXAMPLE OF USE (copy and paste the code in the command line) 
%   Random data, two factors, with 4 and 3 levels, but only the first one 
%   is significative, and 4 replicates
%
% reps = 4;
% vars = 400;
% levels = {[1,2,3,4],[1,2,3]};
% 
% F = create_design(levels,'Replicates',reps);
% 
% X = zeros(size(F,1),vars);
% for ii = 1:length(levels{1}),
%     X(find(F(:,1) == levels{1}(ii)),:) = simuleMV(length(find(F(:,1) == levels{1}(ii))),vars,'LevelCorr',8) + repmat(randn(1,vars),length(find(F(:,1) == levels{1}(ii))),1);
% end
% 
% table = parglm(X, F)
%
%
% EXAMPLE OF USE (copy and paste the code in the command line)
%   Random data, two significative factors, with 4 and 3 levels, and 4 
%   replicates
%
% reps = 4;
% vars = 400;
% levels = {[1,2,3,4],[1,2,3]};
% 
% F = create_design(levels,'Replicates',reps);
% 
% X = zeros(size(F,1),vars);
% for ii = 1:length(levels{1}),
%     fi{ii} = randn(1,vars);
% end
% for jj = 1:length(levels{2}),
%     fj{jj} = randn(1,vars);
% end
% for ii = 1:length(levels{1}),
%     for jj = 1:length(levels{2}),
%         X(find(F(:,1) == levels{1}(ii) & F(:,2) == levels{2}(jj)),:) = simuleMV(reps,vars,8) + repmat(fi{ii} + fj{jj},reps,1);
%     end
% end
% 
% table = parglm(X, F, 'Model',{[1 2]})
%
%
% EXAMPLE OF USE (copy and paste the code in the command line)
%   Random data, two factors with 4 and 3 levels, and 4 replicates, with 
%   significant interaction
%
% reps = 4;
% vars = 400;
% levels = {[1,2,3,4],[1,2,3]};
%
% F = create_design(levels,'Replicates',reps);
%
% X = zeros(size(F,1),vars);
% for ii = 1:length(levels{1}),
%     for jj = 1:length(levels{2}),
%         X(find(F(:,1) == levels{1}(ii) & F(:,2) == levels{2}(jj)),:) = simuleMV(reps,vars,8) + repmat(randn(1,vars),reps,1);
%     end
% end
%
% table = parglm(X, F, 'Model',{[1 2]})
%
%
% coded by: Jos� Camacho (josecamacho@ugr.es)
% last modification: 23/Apr/2024
%
% Copyright (C) 2024  Universidad de Granada
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

%% Arguments checking

% Set default values
routine=dbstack;
assert (nargin >= 2, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);
N = size(X, 1);
M = size(X, 2);

n_factors = size(F,2);                 % number of factors

% Introduce optional inputs as parameters (name-value pair) 
p = inputParser;

addParameter(p,'Model','linear'); 
addParameter(p,'Preprocessing',2);
addParameter(p,'Permutations',1000);  
addParameter(p,'Ts',1); 
addParameter(p,'Ordinal',zeros(1,size(F,2))); 
addParameter(p,'Fmtc',0); 
addParameter(p,'Coding',zeros(1,size(F,2))); 
addParameter(p,'Nested',[]); 
parse(p,varargin{:});

% Extract inputs from inputParser for code legibility
model = p.Results.Model;
prep = p.Results.Preprocessing;
n_perm = p.Results.Permutations;
ts = p.Results.Ts;
ordinal = p.Results.Ordinal;
fmtc = p.Results.Fmtc;
coding = p.Results.Coding;
nested = p.Results.Nested;

if isempty(model), model = 'linear'; end;

if isequal(model,'linear')
    interactions = [];
end  
    
f = 1:n_factors;
if ~isempty(nested), f(nested(:,2)) = []; end
if isequal(model,'interaction')
    interactions = allinter(f,2);
end    

if isequal(model,'full')
    interactions = allinter(f,length(f));
end    

if isnumeric(model) && isscalar(model) && model >= 2 && model <= n_factors
        interactions = allinter(f,model);
end    

if isnumeric(model) && ~isscalar(model)
    for ii = 1:size(model,1)
        interactions{ii} = model(ii,:);
    end
end    

if iscell(model), interactions = model; end 

% Validate dimensions of input data
assert (isequal(size(prep), [1 1]), 'Dimension Error: parameter ''Preprocessing'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(n_perm), [1 1]), 'Dimension Error: parameter ''Permutations'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(ts), [1 1]), 'Dimension Error: parameter ''Ts'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(ordinal), [1 size(F,2)]), 'Dimension Error: parameter ''Ordinal'' must be 1-by-F. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(fmtc), [1 1]), 'Dimension Error:parameter ''Fmtc'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(coding), [1 size(F,2)]), 'Dimension Error: parameter ''Coding'' must be 1-by-F. Type ''help %s'' for more info.', routine(1).name);


%% Main code
                  
n_interactions      = length(interactions);      % number of interactions
n_factors           = size(F,2);                 % number of factors
if fmtc
    mtcc                = n_factors + n_interactions;        % correction for the number of tests
else
    mtcc = 1;
end
SSQ_factors         = zeros(n_perm*mtcc,n_factors,1);      % sum of squares for factors
SSQ_interactions    = zeros(n_perm*mtcc,n_interactions);   % sum of squares for interactions
F_factors           = zeros(n_perm*mtcc + 1,n_factors,1);      % F for factors
F_interactions      = zeros(n_perm*mtcc + 1,n_interactions);   % F for interactions
p_factor            = zeros(1, n_factors);       % p-values factors
p_interaction       = zeros(1, n_interactions);  % p-values interactions

% check if data is complex
iscomp              = ~isreal(X);

% In column space
parglmo.factors                = cell(n_factors,1);
parglmo.interactions           = cell(n_interactions,1);

% preprocess the data
[Xs,m,dt] = preprocess2D(X,'Preprocessing',prep);
X = X./(ones(size(X,1),1)*dt);

% Make structure with general 'variables'
parglmo.data            = X;
parglmo.prep            = prep;
parglmo.scale           = dt;
parglmo.design          = F;
parglmo.n_factors       = n_factors;
parglmo.n_interactions  = n_interactions;
parglmo.n_perm          = n_perm;
parglmo.ts              = ts;
parglmo.ordinal         = ordinal;
parglmo.fmtc            = fmtc;
parglmo.coding          = coding;
parglmo.nested          = nested; 

% Create Design Matrix
if prep
    n = 1;
    D = ones(size(X,1),1);
else
    n = 0;
    D = [];
end

for f = 1 : n_factors
    if ordinal(f)
        D(:,n+1) = preprocess2D(F(:,f),'Preprocessing',1);
        parglmo.factors{f}.Dvars = n+1;
        n = n + 1;
        parglmo.factors{f}.order = 1;
        parglmo.factors{f}.factors = [];
    else
        if isempty(nested) || isempty(find(nested(:,2)==f)) % if not nested
            parglmo.factors{f}.factors = [];
            uF = unique(F(:,f));
            parglmo.n_levels(f) = length(uF);
            for ii = 2:length(uF)
                D(find(ismember(F(:,f),uF(ii))),n+ii-1) = 1;
            end
            parglmo.factors{f}.Dvars = n+(1:length(uF)-1);
            if coding(f) == 1
                D(find(ismember(F(:,f),uF(1))),parglmo.factors{f}.Dvars) = 0;
            else
                D(find(ismember(F(:,f),uF(1))),parglmo.factors{f}.Dvars) = -1;
            end
            n = n + length(uF) - 1;
            parglmo.factors{f}.order = 1;
        else % if nested
            ind = find(nested(:,2)==f);
            ref = nested(ind,1);
            parglmo.factors{f}.factors = [ref parglmo.factors{ref}.factors];
            urF = unique(F(:,ref));
            parglmo.n_levels(f) = 0;
            parglmo.factors{f}.Dvars = [];
            for jj = 1:length(urF)
                rind = find(ismember(F(:,ref),urF(jj)));
                uF = unique(F(rind,f));
                parglmo.n_levels(f) = parglmo.n_levels(f) + length(uF);
                for ii = 2:length(uF)
                    D(rind(find(ismember(F(rind,f),uF(ii)))),n+ii-1) = 1;
                end
                parglmo.factors{f}.Dvars = [parglmo.factors{f}.Dvars n+(1:length(uF)-1)];
                if coding(f) == 1
                    D(rind(find(ismember(F(rind,f),uF(1)))),n+(1:length(uF)-1)) = 0;
                else
                    D(rind(find(ismember(F(rind,f),uF(1)))),n+(1:length(uF)-1)) = -1;
                end
                n = n + length(uF) - 1;
            end   
            parglmo.factors{f}.order = parglmo.factors{ref}.order + 1;
        end
    end
end

for ii = 1 : n_interactions
    Dout = computaDint(interactions{ii},parglmo.factors,D);
    D = [D Dout];
    parglmo.interactions{ii}.Dvars = n+1:size(D,2);
    parglmo.interactions{ii}.factors = interactions{ii};
    n = size(D,2);
    parglmo.interactions{ii}.order = max(parglmo.factors{interactions{ii}(1)}.order,parglmo.factors{interactions{ii}(2)}.order) + 1;
end

% Degrees of freedom
Tdf = size(X,1); 
if prep
    mdf = 1;
else
    mdf = 0;
end
Rdf = Tdf-mdf;
for f = 1 : n_factors
    if ordinal(f)
        df(f) = 1;
    else
        df(f) = length(parglmo.factors{f}.Dvars);
    end
    Rdf = Rdf-df(f);
end
df_int = [];
for ii = 1 : n_interactions
    df_int(ii) = prod(df(parglmo.interactions{ii}.factors));
    Rdf = Rdf-df_int(ii);
end
if Rdf < 0
    disp('Warning: degrees of freedom exhausted');
    return
end

% Handle missing data 
[r,c]=find(isnan(X));
Xnan = X;
ru = unique(r);
for ii=1:length(ru)
    ind = find(r==ru(ii));
    ind2 = find(sum((D-ones(size(D,1),1)*D(r(ind(1)),:)).^2,2)==0);
    for jj=1:length(c(ind))
        ind3 = find(isnan(X(ind2,c(ind(jj)))));
        if length(ind2)>length(ind3)
            X(r(ind(jj)),c(ind(jj))) = nanmean(X(ind2,c(ind(jj)))); % use conditional mean replacement
        else
            X(r(ind(jj)),c(ind(jj))) = nanmean(X(:,c(ind(jj)))); % use unconditional mean replacement if CMR not possible
        end
    end
end
parglmo.data = X;
parglmo.Xnan = Xnan;

if iscomp
    SSQ_X = sum(abs(diag(X*X.')));
else
    SSQ_X = sum(sum((X).^2));
end

% GLM model calibration with LS, only fixed factors
pD =  pinv(D'*D)*D';
B = pD*X;
X_residuals = X - D*B;
parglmo.D = D;
parglmo.B = B;
%parglmo.mean = parglmo.D*parglmo.B(:,1);

% Create Effect Matrices
if prep
    parglmo.inter = D(:,1)*B(1,:);
    if iscomp
        SSQ_inter = sum(abs(diag(parglmo.inter * parglmo.inter.')));
    else
        SSQ_inter = sum(sum((parglmo.inter).^2));
    end
else
    parglmo.inter = 0;
    SSQ_inter = 0;
end    
%SSQ_residuals = sum(sum((X_residuals).^2));

if iscomp
    SSQ_residuals = sum(abs(diag(X_residuals*X_residuals.')));
else
    SSQ_residuals = sum(sum((X_residuals).^2));
end

for f = 1 : n_factors
    parglmo.factors{f}.matrix = D(:,parglmo.factors{f}.Dvars)*B(parglmo.factors{f}.Dvars,:);
    if iscomp
        SSQ_factors(1,f) = sum(abs(diag(parglmo.factors{f}.matrix * parglmo.factors{f}.matrix.')));
    else
        SSQ_factors(1,f) = sum(sum((parglmo.factors{f}.matrix).^2)); % Note: we are not using Type III sum of squares, and probably we should, although we did not find any difference in our experiments
    end
end

for ii = 1 : n_interactions
    parglmo.interactions{ii}.matrix = D(:,parglmo.interactions{ii}.Dvars)*B(parglmo.interactions{ii}.Dvars,:);
    if iscomp
        SSQ_interactions(1,ii) = sum(abs(diag(parglmo.interactions{ii}.matrix * parglmo.interactions{ii}.matrix.')));
    else
        SSQ_interactions(1,ii) = sum(sum((parglmo.interactions{ii}.matrix).^2));
    end
end

for f = 1 : n_factors
    if ts==2 % pooled reference MSE
        SS_ref = 0;
        Df_ref = 0;
        for f2 = 1 : n_factors
            if ~isempty(find(f==parglmo.factors{f2}.factors))
                SS_ref = SS_ref + SSQ_factors(1,f2);
                Df_ref = Df_ref + df(f2);
            end
        end
        for ii = 1 : n_interactions
            if ~isempty(find(f==parglmo.interactions{ii}.factors))
                SS_ref = SS_ref + SSQ_interactions(1,ii);
                Df_ref = Df_ref + df_int(ii);
            end
        end
        if SS_ref == 0
            F_factors(1,f) = (SSQ_factors(1,f)/df(f))/(SSQ_residuals/Rdf);
        else
            F_factors(1,f) = (SSQ_factors(1,f)/df(f))/(SS_ref/Df_ref);
        end
    else
        F_factors(1,f) = (SSQ_factors(1,f)/df(f))/(SSQ_residuals/Rdf);
    end
end

for ii = 1 : n_interactions
    F_interactions(1,ii) = (SSQ_interactions(1,ii)/df_int(ii))/(SSQ_residuals/Rdf);
end

parglmo.effects = 100*([SSQ_inter SSQ_factors(1,:) SSQ_interactions(1,:) SSQ_residuals]./SSQ_X);
parglmo.residuals = X_residuals;

% Permutations
parfor jj = 1 : n_perm*mtcc
   
    perms = randperm(size(Xnan,1)); % permuted data (permute whole data matrix)
    
    X = Xnan(perms, :);
    [r,c]=find(isnan(X));
    ru = unique(r);
    for ii=1:length(ru)
        ind = find(r==ru(ii));
        ind2 = find(sum((D-ones(size(D,1),1)*D(r(ind(1)),:)).^2,2)==0);
        for f=1:length(c(ind))
            ind3 = find(isnan(X(ind2,c(ind(f)))));
            if length(ind2)>length(ind3)
                X(r(ind(f)),c(ind(f))) = nanmean(X(ind2,c(ind(f)))); % use conditional mean replacement
            else
                X(r(ind(f)),c(ind(f))) = nanmean(X(:,c(ind(f)))); % use unconditional mean replacement if CMR not possible
            end
        end
    end
      
    B = pD*X;
    X_residuals = X - D*B;
    if iscomp
        SSQ_residualsp = sum(abs(diag(X_residuals*X_residuals.')));
    else
        SSQ_residualsp = sum(sum((X_residuals).^2));
    end

    % Factors
    factors = cell(1,n_factors);
    SSQf = zeros(1,n_factors);
    for f = 1 : n_factors
        factors{f}.matrix = D(:,parglmo.factors{f}.Dvars)*B(parglmo.factors{f}.Dvars,:);
        if iscomp
            SSQf(f) = sum(abs(diag(factors{f}.matrix*factors{f}.matrix.')));
        else
            SSQf(f) = sum(sum((factors{f}.matrix).^2));
        end
    end
    SSQ_factors(1 + jj,:) = SSQf;
    
    % Interactions
    interacts = {};
    SSQi = zeros(1,n_interactions);
    for ii = 1 : n_interactions
        interacts{ii}.matrix = D(:,parglmo.interactions{ii}.Dvars)*B(parglmo.interactions{ii}.Dvars,:);    
        if iscomp
            SSQi(ii) = sum(abs(diag(interacts{ii}.matrix * interacts{ii}.matrix.')))
        else
            SSQi(ii) = sum(sum((interacts{ii}.matrix).^2));
        end
    end
    SSQ_interactions(1 + jj,:) = SSQi;
    
    % F Factors
    Ff = zeros(1,n_factors);
    for f = 1 : n_factors
        if ts==2% pooled reference MSE
            MSS_ref = 0;
            Df_ref = 0;
            for f2 = 1 : n_factors
                if ~isempty(find(f==parglmo.factors{f2}.factors))
                    MSS_ref = MSS_ref + SSQf(f2);
                    Df_ref = Df_ref + df(f2);
                end
            end
            for ii = 1 : n_interactions
                if ~isempty(find(f==parglmo.interactions{ii}.factors))
                    MSS_ref = MSS_ref + SSQi(ii);
                    Df_ref = Df_ref + df_int(ii);
                end
            end
            if MSS_ref == 0
                Ff(f) = (SSQf(f)/df(f))/(SSQ_residualsp/Rdf);
            else
                Ff(f) = (SSQf(f)/df(f))/(MSS_ref/Df_ref);
            end
        else
            Ff(f) = (SSQf(f)/df(f))/(SSQ_residualsp/Rdf);
        end
    end
    F_factors(1 + jj,:) = Ff; 
    
    % F Interactions
    Fi = zeros(1,n_interactions);
    for ii = 1 : n_interactions
        Fi(ii) = (SSQi(ii)/df_int(ii))/(SSQ_residualsp/Rdf);
    end
    F_interactions(1 + jj,:) = Fi; 

end        

% Select test statistic
if ts
    ts_factors = F_factors;
    ts_interactions = F_interactions;
else
    ts_factors = SSQ_factors;
    ts_interactions = SSQ_interactions;
end
    
% Calculate p-values
for f = 1 : n_factors
    p_factor(f) = (size(find(ts_factors(2:(n_perm*mtcc + 1), f) >= ts_factors(1, f)),1) + 1)/(n_perm*mtcc+1);
end
for ii = 1 : n_interactions
    p_interaction(ii) = (size(find(ts_interactions(2:(n_perm*mtcc + 1), ii) ...
        >= ts_interactions(1, ii)),1) + 1)/(n_perm*mtcc+1);
end

% Multiple test correction for several factors/interactions
parglmo.p = [p_factor p_interaction]; 
if mtcc > 1
    switch fmtc
        case 1 % Bonferroni 
            parglmo.p = min(1,parglmo.p * mtcc); 

        case 2 % Holm/Hochberg
            [~,indx] = sort(parglmo.p,'ascend');
            for ind = 1 : mtcc 
                parglmo.p(indx(ind)) = min(1,parglmo.p(indx(ind)) * (mtcc-ind+1));
            end

        case 3 % Benjamini & Hochberg
            [~,indx] = sort(parglmo.p,'ascend');
            parglmo.p(indx(mtcc)) = parglmo.p(indx(mtcc));
            for ind = mtcc-1 : -1 : 1 
                parglmo.p(indx(ind)) = min([1,parglmo.p(indx(ind)) * mtcc/ind]);
            end

        case 4 % Q-value from Benjamini & Hochberg
            [~,indx] = sort(parglmo.p,'ascend');
            parglmo.p(indx(mtcc)) = parglmo.p(indx(mtcc));
            for ind = mtcc-1 : -1 : 1 
                parglmo.p(indx(ind)) = min([1,parglmo.p(indx(ind)) * mtcc/ind,parglmo.p(indx(ind+1))]);
            end
    end
end


%% ANOVA-like output table

name={'Mean'};
for f = 1 : n_factors
    name{end+1} = sprintf('Factor %d',f);
end
for ii = 1 : n_interactions
    name{end+1} = sprintf('Interaction %s',strrep(num2str(parglmo.interactions{ii}.factors),'  ','-'));
end
name{end+1} = 'Residuals';
name{end+1} = 'Total';
   

SSQ = [SSQ_inter SSQ_factors(1,:) SSQ_interactions(1,:) SSQ_residuals SSQ_X];
par = [parglmo.effects 100];
DoF = [mdf df df_int Rdf Tdf];
MSQ = SSQ./DoF;
F = [nan F_factors(1,:) F_interactions(1,:) nan nan];
p_value = [nan parglmo.p nan nan];

%T = table(name', SSQ', par', DoF', MSQ', F', p_value','VariableNames', {'Source','SumSq','PercSumSq','df','MeanSq','F','Pvalue'});
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
if isOctave
  T.mat = [SSQ', par', DoF', MSQ', F', p_value'];
  T.var = {'SumSq', 'PercSumSq', 'df', 'MeanSq', 'F', 'Pvalue'};
  T.source = name';
else
  T = table(name', SSQ', par', DoF', MSQ', F', p_value','VariableNames', {'Source','SumSq','PercSumSq','df','MeanSq','F','Pvalue'});
end

end

%% Auxiliary function for interactions

function interactions = allinter(factors,order)
    
    if order > 2
        interactions = allinter(factors,order-1);
        for ii = 1:length(interactions)
            for jj = factors(find(factors > max(interactions{ii})))
                interactions{end+1} = [interactions{ii} jj];
            end
        end
    else
        interactions = {};
        for ii = factors
            for jj = factors(find(factors >ii))
                interactions{end+1} = [ii jj];
            end
        end
    end
    
end
    
        
function Dout = computaDint(interactions,factors,D) % Compute coding matrix

    if length(interactions)>1
        deepD = computaDint(interactions(2:end),factors,D);
        Dout = [];
        for k = factors{interactions(1)}.Dvars
            for l = 1:size(deepD,2)
                Dout(:,end+1) = D(:,k).* deepD(:,l);
            end
        end
    else
        Dout = D(:,factors{interactions}.Dvars);
    end

end

