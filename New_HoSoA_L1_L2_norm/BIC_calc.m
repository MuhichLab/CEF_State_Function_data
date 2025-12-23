clear;
clc;
%% Load errors


model = "BSF_model";

load(model +'_Observed_dds.mat')
load(model +'_Predicted_dds.mat')
load(model +'_Predicted_MU.mat')
load(model +'_Observed_MU.mat')
% load(model,'numLs');

numLs = 22

%% Bayesian information criterion paramter setup

numParam = transpose(3:3:numLs*3)+6;

%% Deterministic Model 
% % i.e. for any specific choice of model parameters there is no uncertainty
% % in the predicted value.
% % NOTE: var is std.^2
% 
%     [~,numObs] = size(observed_exp); % number of total observations
%     obsv = observed_exp(1,:);
%     sig_errs = std((predicted_exp-obsv),[],2);
%     model_error = ((predicted_exp-obsv)).^2;
%     sum_error = sum(model_error,2);
% 
%     logLs = -numObs/2*log(2*pi()) - numObs*log(sig_errs) - (1./(2*sig_errs.^2)).*sum_error;


    %% Deterministic Model 
% i.e. for any specific choice of model parameters there is no uncertainty
% in the predicted value.
% NOTE: var is std.^2
% 
    % [~,numObs] = size(observed_exp); % number of total observations
    % obsv = observed_exp(1,:);
    % % sig_errs = std((predicted_exp-obsv),[],2);
    % model_errors = (predicted_exp-obsv);
    % mean_errors = mean(model_errors,2);
    % sum_errors = sum((model_errors - mean_errors).^2,2);
    % sig_errs = std((model_errors - mean_errors),[],2);
    % 
    % logLs = -numObs/2*log(2*pi()) - numObs*log(sig_errs) - (1./(2*sig_errs.^2)).*sum_errors;

%% Probabilistic Model
%     % i.e. The model prediction includes a statistical noise componenet (the
%     % noise is in the data)
% 
numObs = 100; % number of total observations
comp_model_err = predicted_exp - observed_exp;
model_err = [];

count = 1;
sigma_data = [];
for s = 22:22:2200
    model_err(:,count) = mean(predicted_exp(:,s-21:s),2) - mean(observed_exp(1,s-21:s));
    sigma_data(1,count) = std(observed_exp(1,s-21:s));
    count = count + 1;
end

sigma_err = std(model_err,[],2);
sq_err = model_err.^2;
sum_sigs = sigma_data.^2 + sigma_err.^2;

% Log Likelihood
logLs = -numObs/2*log(2*pi()) - 1/2*sum(log(sum_sigs),2) - 1/2*sum(sq_err/sum_sigs,2);
    
%% BIC Calc using Matlab
[AIC,BIC,ic] = aicbic(logLs,numParam,numObs,Normalize=true);

figure
plot((1:numLs),BIC,'-ok','linewidth',2)
ylabel('BIC')
xlabel('Number of L Excess Terms')

[val,idx] = min(BIC);

best_model = idx

%% 1 to 1 plot


for model = 1:4
    MAE_d = mean(abs(observed_exp(model,:) - predicted_exp(model,:)));
    figure
    hold on
    scatter(observed_exp(model,:),predicted_exp(model,:))
    plot([0 0.35],[0 0.35],':k')
    txt="MAE: " + num2str(sprintf('%.2e',MAE_d));
    text(0.1,0.3,txt)
    xlabel('Observed \delta')
    ylabel('Predicted \delta')
    title("Model " + num2str(model))
    box on
    
    
    MAE_mu = mean(abs(-observed_mu(model,:) - predicted_mu(model,:)));
    figure
    hold on
    scatter(-observed_mu(model,:),predicted_mu(model,:))
    plot([1 2.5],[1 2.5],':k')
    txt="MAE: " + num2str(sprintf('%.2e',MAE_mu));
    text(1.1,2,txt)
    xlabel('Observed \mu')
    ylabel('Predicted \mu')
    title("Model " + num2str(model))
    box on
    
end

