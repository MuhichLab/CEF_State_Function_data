clear;
clc;

%% Build CEF

syms x y z T To real
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



dft_data = table2array(readtable('../BSF_Full_DFT_set.xlsx'));
X =  dft_data(:,1); % mol frac x
Y =  dft_data(:,2); % mol frac y
Z =  dft_data(:,3); % delta
eV = dft_data(:,4); % E eV per ABO3

% Set the H(T=0) values
% gs order is: (must be in this order)
%(1) ABO3
%(2) A'BO3
%(3) AB'O3
%(4) A'B'O3
%(5) ABO2.5
%(6) A'BO2.5
%(7) AB'O2.5
%(8) A'B'O2.5
gs = zeros(1,8);
gs(1) = eV(X==0&Y==0&Z==0);
gs(2) = eV(X==1&Y==0&Z==0);
gs(3) = eV(X==0&Y==1&Z==0);
gs(4) = eV(X==1&Y==1&Z==0);
gs(5) = eV(X==0&Y==0&Z==0.5);
gs(6) = eV(X==1&Y==0&Z==0.5);
gs(7) = eV(X==0&Y==1&Z==0.5);
gs(8) = eV(X==1&Y==1&Z==0.5);

% the chemical eqution form we are testing is: 
% A(1-x)A'xBO3-z
xc = 1; % Means there is a substiution on the A lattice
yc = 0; % Means there is no substiution on the B lattice

% Cp expanion is constant ( cp = a )
cp = 1; % Means one term in the heat cpacity expansion

% Create General CEF Function
[G_soln,G_ex,G_o] = CEF(xc,yc);

% Simplify excess terms - remove parallel terms
[G_soln,G_ex,subin,subout] = Simplify_Excess(G_soln,G_ex,xc,yc);

% Expand Gibbs endmember terms to the for G(T) = A + B*T + C*T*ln(T) ...
[G_sol,G_dft,Go,Go_dft] = T_Expand_endmembers_w_diff(G_soln,G_o,cp);

% Expand Gibbs excess terms to the for G(T) = A + B*T + C*T*ln(T) ...
[Gsol,Gdft,Gex] = T_Expand_excess(G_sol,G_dft,G_ex,cp);

% Since we have DFT data we can "localize" the the enthalpic energy
[Gsol,Go,Gdft,Go_dft] = Localize_H_w_diff(Gsol,Go,Gdft,Go_dft,gs);

% Make term labeleing sequential
[Gsol, Gdft,ex_terms, end_terms, ex_terms_og, end_terms_og, numLs]...
    = Get_Seq_Terms_w_diff(Gsol,G_soln,Gdft,xc,yc);

save("Model_Workspace");