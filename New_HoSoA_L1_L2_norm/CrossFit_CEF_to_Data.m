clear;
clc;

% set seed for reproducibility
rng(42)
start  = tic();


syms x y z T real
syms G [1 8] real
syms gDiffA [1 8] real
syms gDiffB [1 8] real
syms gDiffC [1 8] real
syms gDiffD [1 8] real
syms gDiffHo [1 8] real
syms gDiffSo [1 8] real
syms Ho [1 80] real
syms So [1 80] real
syms A [1 80] real
syms B [1 80] real
syms C [1 80] real
syms D [1 80] real
syms L [1 80] real

%% Load CEF model wroskpace
load("CEF_Model_generator/Model_Workspace.mat")
num_ex_params = length(ex_terms_og);
num_end_params = length(end_terms_og);

dG_dy = diff(Gsol,z);

%% some job settings
file_save_name = 'fminunc_1e-2';
val = 1e-2;  % intial guess value

% how much do we weight the DFT data?? 
% first number is DFT weight second is exp weight
lambda = [1,1]; %using ratio of number of data points

% L1 L2 norm Regularization coefficients
lambda_L1 = 1e-2; % Coefficient for L1 norm # Higher values increase sparsity in the parameter set.
lambda_L2 = 1e-3; % Coefficient for L2 norm # Higher values enforce small parameter magnitudes.

%% Load Data

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

% x/y_vals are the unique experimental mol fractions gathered
x_vals = unique(x_exp);
y_vals = unique(y_exp);

    % X0 = is the instial guess for the d0 solvers need for first fit
X0 = zeros(length(x_vals),1)+0.05; 
Y0 = zeros(length(x_vals),1)+0.05;

%% Calculate O chemical potential at TGA data points

H_o = zeros(length(Temp),1);
S_o = zeros(length(Temp),1);
for t = 1:length(Temp)
    [H_o2,S_o2] = get_O2_thermo(Temp(t));
    H_o(t) = H_o2;
    S_o(t) = S_o2;
end
conv = 96.487;  % 1 eV = 96.487 kJ/mol
R_gas = 8.314462/1000/conv; %eV/K

muhg_o2 = H_o/conv - (Temp).*S_o/conv + R_gas.*(Temp).*log(Press);

muhg_o = 0.5*muhg_o2;
%% If delta values are a delta-delta based on a refrence point set that here

T_ref = 573.15;  % K
P_ref = 0.9; % Bar
d_chem = round(get_mu_o(T_ref,P_ref),4);

% % Set to 0 if delta values are aboslute
% T_ref = 0;  % K
% P_ref = 0;  % Bar
% d_chem = 0;

dref = zeros(1,length(x_vals)*length(y_vals));
%%
ncores = feature('numcores');
paracores = ncores/2;

poolobj = gcp('nocreate'); % If no pool, do not create new one.
if isempty(poolobj)
    poolobj = parpool('local',paracores);
end
   
%% Set the intial guess values sor the parameter optimizer

% Ao Ho So are hard coded need to fix
add_me_end = zeros(1,num_end_params);
add_me_ex = zeros(1,num_ex_params);

A0s_end = ismember(end_terms,gDiffA);
A0s_ex = ismember(ex_terms,A);

H0s_end = ismember(end_terms,gDiffHo);
H0s_ex = ismember(ex_terms,Ho);

S0s_end = ismember(end_terms,gDiffSo);
S0s_ex = ismember(ex_terms,So);

add_me_end(A0s_end) = val;
add_me_end(H0s_end) = val;
add_me_end(S0s_end) = val;

%add some noise
add_me_end = add_me_end + rand(size(add_me_end))/(1/(val/10));

% add_me_end(B0s_end) = 1e-4;

add_me_ex(A0s_ex) = val;
add_me_ex(H0s_ex) = val;
add_me_ex(S0s_ex) = val;

%add some noise
add_me_ex = add_me_ex + rand(size(add_me_ex))/(1/(val/10));

%% fmincon options
% Uncomment the 'PlotFcn if you wish to watch the minimization while
% running on your local machine - BE CAREFUL of parentheses

options = optimoptions('fminunc','MaxFunctionEvaluations',1000000,...
    'MaxIterations',1000000,'StepTolerance',1e-9,'UseParallel',true,...
    'Display','off','OptimalityTolerance',1e-6,'FunctionTolerance',1e-6,'PlotFcn',{'optimplotfval'});

%% Endmebers only fit - i.e Go parameter minimization 

s1 = tic();
% Endmember self intilization
guess_excess = zeros(1,num_ex_params);
guess0 = zeros(1,num_end_params)+add_me_end;
dp_terms = ones(num_ex_params,1);


% the intial guesss with an unknown dref is really really bad. First we
% optiize the paramters off a dref guess that does not change. Then they
% are reoptmized in a self consistent loop with dref.
if T_ref ~= 0
    % NOTE: I pass in -muhg_o because the fit is dG/dy = -mu
    [guess_end,fval,~,output1] = ...
        fminunc(@(guess_end)cross_obj_fun_sym(guess_end,guess_excess,X,Y,Z,eV,...
        x_exp,y_exp,dd,X0,Y0,x_vals,y_vals,-muhg_o,Temp,Gdft,dG_dy,d_chem,...
        lambda,T_ref,dp_terms,num_ex_params,false,ex_terms,end_terms,lambda_L1,lambda_L2,...
        x,y,z,T,gDiffA,gDiffB,gDiffC,gDiffD,gDiffHo,gDiffSo,Ho,So,A,B,C,D)...
        ,guess0,options); 
    
    guess0 = guess_end;
    fprintf('\nThe # of func calls for the Endmembers with fixed do was : %0.5f \n',output1.funcCount);

end

[guess_end,end_error,~,output2] = ...
    fminunc(@(guess_end)cross_obj_fun_sym(guess_end,guess_excess,X,Y,Z,eV,...
    x_exp,y_exp,dd,X0,Y0,x_vals,y_vals,-muhg_o,Temp,Gdft,dG_dy,d_chem,...
    lambda,T_ref,dp_terms,num_ex_params,true,ex_terms,end_terms,lambda_L1,lambda_L2,...
    x,y,z,T,gDiffA,gDiffB,gDiffC,gDiffD,gDiffHo,gDiffSo,Ho,So,A,B,C,D)...
    ,guess0,options); 

fprintf('\nThe RSS for the Endmembers only fit: %0.5f \n',end_error);
fprintf('\nThe # of func calls for the Endmembers with solve do was : %0.5f \n',output2.funcCount);


if T_ref ~= 0
    dGdy = subs(dG_dy,ex_terms,guess_excess);
    dGdy = matlabFunction(subs(dGdy,end_terms,guess_end));
    d0_end_only = solve_dref(dGdy,x_vals,d_chem,T_ref);
else
    d0_end_only = zeros(length(x_vals),1);
end


%% Endmember values for in progress watching

final_tab_endMem = array2table(guess_end', ...
      'VariableNames',{'Coeff'});
final_tab_endMem.Properties.RowNames = string(end_terms');

%% Logic dropping of terms and refitting

% initilize excess terms
guess_hold = zeros(1,num_ex_params)+add_me_ex;
% guess_excess = zeros(1,num_ex_params)+add_me_ex;
%     guess_excess = fake_excess_terms;
% initilize some populated arrays
dp_terms = ones(1,length(guess_excess)); % dp_terms controls which terms to set to 0 effectivily dropping them
d0_progress = zeros(length(x_vals),numLs); % tracks dref calculations
all_fits = zeros(length(guess_excess)+length(guess_end),numLs); % tracks optimized parameters
spot_track = zeros(1,numLs); % used for tracking of other arrays
count = 21;
num_terms = numLs+1; % arbitrary value higher than the number of excess terms
RSS = zeros(1,numLs);

% Fit L17 and L18 only
dp_terms(1:16) = 0;
dp_terms(1+22:16+22) = 0;
dp_terms(1+44:16+44) = 0;
dp_terms(19:22) = 0;
dp_terms(19+22:22+22) = 0;
dp_terms(19+44:22+44) = 0;

%%%%%% Either reset paramaters or use parameters from last ittr %%%%%%%%%%%
% %     keeps previous values and removes dropped terms
%    guess0 = guess_excess(dp_terms>0);
%     % reset terms to intital guess
    guess0 = guess_hold(dp_terms>0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    optT1 = tic();
    
    % Again the intial guesss with an unknown dref is really really bad.
    % This time we can use the endmember dref already optimized for
    if T_ref ~= 0  
        % NOTE: I pass in -muhg_o because the fit is dG/dy = -mu
        [guess_excess,~,~,output3]  =...
            fminunc(@(guess_excess)cross_obj_fun_sym(guess_end,guess_excess,X,Y,Z,eV,...
        x_exp,y_exp,dd,d0_end_only,Y0,x_vals,y_vals,-muhg_o,Temp,Gdft,dG_dy,d_chem,...
        lambda,0,dp_terms,num_ex_params,false,ex_terms,end_terms,lambda_L1,lambda_L2,...
        x,y,z,T,gDiffA,gDiffB,gDiffC,gDiffD,gDiffHo,gDiffSo,Ho,So,A,B,C,D),guess0,options);

        fprintf('\nThe # of func calls for the excess with fixed do was : %0.5f \n',output3.funcCount);

    end
    
    guess0 = guess_excess;
    
    [guess_excess,drop_error,exitflag,output4]  =...
        fminunc(@(guess_excess)cross_obj_fun_sym(guess_end,guess_excess,X,Y,Z,eV,...
    x_exp,y_exp,dd,d0_end_only,Y0,x_vals,y_vals,-muhg_o,Temp,Gdft,dG_dy,d_chem,...
    lambda,T_ref,dp_terms,num_ex_params,true,ex_terms,end_terms,lambda_L1,lambda_L2,...
    x,y,z,T,gDiffA,gDiffB,gDiffC,gDiffD,gDiffHo,gDiffSo,Ho,So,A,B,C,D),guess0,options);
    fprintf('\nOptimization took (mins): %0.2f \n',toc(optT1)/60);
    fprintf('\nThe final RSS error for this itteration: %0.5f \n',drop_error);
    RSS(count) = drop_error;
    fprintf('\nThe # of func calls for the excess with solve do was : %0.5f \n',output4.funcCount);

    % rebuild the full guess_drop with zeros in place for dropped terms - THIS IS necessary here!
    og_guess = guess_excess;
    guess_excess = zeros(1,length(dp_terms));
    n = 1;
    for k = (1:length(dp_terms))
        if dp_terms(k) == 1
            guess_excess(k) = og_guess(n);
            n = n + 1;
        end
    end  

    all_fits(:,count) = [guess_excess'; guess_end'];

    if T_ref ~= 0
        dGdy = subs(dG_dy,ex_terms,guess_excess);
        dGdy = matlabFunction(subs(dGdy,end_terms,guess_end));
        d0 = solve_dref(dGdy,x_vals,d_chem,T_ref);
    else
        d0 = zeros(length(x_vals),1);
    end

    d0_progress(:,count) = d0; 

total_time = toc(start);
fprintf('\nTime test for optimzation of all paramters: %0.2f \n',total_time/60);


% %% Table of all fits
% 
% all_tab = array2table(all_fits);
% 
% all_tab.Properties.RowNames = string([ex_terms';end_terms'])

%% Saving 

save(file_save_name + "_workspace");
save(file_save_name + "_allfits",'all_fits');
save(file_save_name + "_d0s",'d0_progress'); 

fit_col = 2;
Gnew = subs(Gsol,[end_terms ex_terms],[all_fits(end-length(end_terms)+1:end,end-(fit_col)+1)' all_fits(1:length(ex_terms),end-(fit_col)+1)']);
Snew = -diff(Gnew,T);
Hnew = Gnew + T*Snew;
save(file_save_name + "_model",'Gnew');

% [Gnew,Hnew,Snew] = plot_model(file_save_name + "_workspace",2);
% save(file_save_name + "_model",'Gnew');
% 
% Gnew_current = Gnew;
% load('fminunc_1e-6_model');
% compare_models(Gnew_current,Gnew);
