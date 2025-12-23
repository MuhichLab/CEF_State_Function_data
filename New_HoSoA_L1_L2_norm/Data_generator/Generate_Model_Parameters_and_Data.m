clear;
clc;

%% Load CEF model 
load("../CEF_Model_generator/Model_Workspace.mat",'ex_terms')
load("../CEF_Model_generator/Model_Workspace.mat",'end_terms')
load("../CEF_Model_generator/Model_Workspace.mat",'numLs')
load("../CEF_Model_generator/Model_Workspace.mat",'Gsol')
load("../CEF_Model_generator/Model_Workspace.mat",'gs')

%% Create Fake Paramter Values From loaded model

[fake_excess_terms] = fake_excess(ex_terms);
[fake_end_terms] = fake_endmembers_w_diff(gs,end_terms,T0);

% random L terms select
combos = nchoosek(1:numLs,3);
ind = randi([1 length(combos)]);
Lterms = combos(ind,:);

fake_excess_terms(setdiff(1:end,[Lterms Lterms+numLs Lterms+(numLs*2)])) = 0;
% fake_excess_terms(1:end) = 0;
%% Create Fake Thermodynamics
Greal = subs(Gsol,[end_terms ex_terms],[fake_end_terms fake_excess_terms]);
Sreal = -diff(Greal,T);
Hreal = Greal + T*Sreal;

thermo_plots(Hreal,Sreal)

%% Create fake Data with noise based on Fake Thermodynamics
% create data pass-in (Model,# of points,noise in delta, # sets of data)
% df is sf with noise
[Tf, pf, xf, df, sf] = create_data(Greal,6,1e-3,2);
Plot_Data(pf,Tf,xf,df)