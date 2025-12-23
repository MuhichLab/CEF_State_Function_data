clear;
clc;

model = "BSF_model";

syms x y z T To real
syms G [1 8] real
syms gDiffA [1 4] real
syms gDiffB [1 4] real
syms gDiffC [1 4] real
syms gDiffD [1 4] real
syms gDiffHo [1 4] real
syms gDiffSo [1 4] real
syms Ho [1 80] real
syms So [1 80] real
syms A [1 80] real
syms B [1 80] real
syms C [1 80] real
syms D [1 80] real
syms L [1 80] real

load(model + "_workspace",'Gsol');
load(model + "_workspace",'Gdft');
load(model + "_workspace",'all_fits');
load(model + "_workspace",'end_terms');
load(model + "_workspace",'ex_terms');
load(model + "_workspace",'numLs');
load(model + "_workspace",'d0_progress');

%% Take necessary derivatives

dG_dy = diff(Gsol,z);

%% Load data
%DFT
dft_data = table2array(readtable('BSF_Full_DFT_set.xlsx'));
rows = any(isnan(dft_data),2);
dft_data(rows,:) = [];
X =  dft_data(:,1); % mol frac x
Y =  dft_data(:,2); % mol frac y
Z =  dft_data(:,3); % delta
eV = dft_data(:,4); % E eV per ABO3

% TGA
exp_data = readtable('Bush_BSF.xlsx');
Press = table2array(exp_data(:,1));  % partial pressure O2 in bar
Temp = table2array(exp_data(:,2));   % K
x_exp = table2array(exp_data(:,3));  %A site fractions
y_exp = table2array(exp_data(:,4));  %B site fractions
dd = table2array(exp_data(:,5));  % delta values

 %% calc dd and mu
num_end_params = length(end_terms);
num_ex_params = length(ex_terms);

% Assing Parameters based on number of excess terms
observed_dft = zeros(numLs,length(eV));
predicted_dft = zeros(numLs,length(eV));
observed_exp = zeros(numLs,length(dd));
predicted_exp = zeros(numLs,length(dd));
observed_mu = zeros(numLs,length(dd));
predicted_mu = zeros(numLs,length(dd));

for spt = 1:numLs
    fit_col = spt

    guess_end = all_fits(end-(num_end_params-1):end,end-(fit_col)+1);
    guess_excess = all_fits(1:num_ex_params,end-(fit_col)+1);

    dGdy = subs(dG_dy,ex_terms,guess_excess');
    dGdy = subs(dGdy,end_terms,guess_end');

%     Gsol = subs(Gsol,ex_terms,guess_excess);
%     Gsol = subs(Gsol,end_terms,guess_end);

    % Gdft = subs(G_dft,intersect(ex_terms,Ho),guess_excess(ismember(ex_terms,Ho))');
    % Gdft = subs(Gdft,intersect(end_terms,gendHo),guess_end(ismember(end_terms,gendHo))');
    % Gdft = subs(Gdft,intersect(ex_terms,Ao),guess_excess(ismember(ex_terms,Ao))');
    % Gdft = subs(Gdft,intersect(end_terms,gendAo),guess_end(ismember(end_terms,gendAo))');

    G_dft = subs(Gdft,intersect(ex_terms,Ho),guess_excess(ismember(ex_terms,Ho))');
    G_dft = subs(G_dft,intersect(end_terms,gDiffHo),guess_end(ismember(end_terms,gDiffHo))');
%     G_dft = subs(Gdft,intersect(ex_terms,A),guess_excess(ismember(ex_terms,A))');
%     G_dft = subs(G_dft,intersect(end_terms,gDiffA),guess_end(ismember(end_terms,gDiffA))');

    parfor n = 1:length(x_exp)

        x_vals = unique(x_exp);
        for k = 1:length(x_vals)
            if x_exp(n) == x_vals(k)
                dref = d0_progress(k,end-(fit_col)+1);
                break
            end
        end
        subbed_dG_soln = matlabFunction(subs(dGdy,[T x],[Temp(n) x_exp(n)]));
        eqn = @(z) (subbed_dG_soln(z)) + get_mu_o(Temp(n),Press(n));
        dd_pred = double(fzero(eqn,[1e-12 0.5-1e-12],optimset('FunValCheck', 'off', 'Display', 'off')));
        %dd_bar = mean(dd(n:n+21));
        observed_exp(fit_col,n) = dd(n);
        predicted_exp(fit_col,n) = dd_pred-dref;
        observed_mu(fit_col,n) = get_mu_o(Temp(n),Press(n));
        predicted_mu(fit_col,n) = subbed_dG_soln(dd(n)+dref);
    end


    for m = 1:length(X)
        subbed_G_soln = matlabFunction(subs(G_dft,T,0));
        observed_dft(fit_col,m) = eV(m);
        predicted_dft(fit_col,m) = double(subbed_G_soln(X(m),Z(m)));
    end
end

save(model +'_Predicted_dds','predicted_exp')
save(model +'_Observed_dds','observed_exp')
save(model +'_Predicted_DFT','predicted_dft')
save(model +'_Observed_DFT','observed_dft')
save(model +'_Predicted_MU','predicted_mu')
save(model +'_Observed_MU','observed_mu')
