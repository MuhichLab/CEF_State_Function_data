function [TT, p, XX, d, s] = create_data(G,N,noise,doop)

    syms x z T muh real
    
    chempot = -diff(G,z) - muh;
    solMu = matlabFunction(chempot);

    T1 = linspace(800,1800,N);
    X = linspace(0,1,N);
    P = logspace(-10,-.6777778,N);
    
    [TT, PP, XX] = meshgrid(T1,P,X);

    TT = reshape(TT,[1 N^3]);
    XX = reshape(XX,[1 N^3]);
    PP = reshape(PP,[1 N^3]);
    
    % make duplicates
    if doop > 1
        pp = PP;
        xx = XX;
        tt = TT;
        for i = 1:doop-1
        PP = [PP,pp];
        XX = [XX,xx];
        TT = [TT,tt];
        end
    end
    
    omuh = solve_for_muO(TT,PP);
    s = 0*omuh;
    for i = 1:length(TT)
        eqn = @(z) solMu(TT(i),omuh(i),XX(i),z);
        s(i) = double(fzero(eqn,[1e-15 0.5-1e-15]));
    end
    
    p = PP;
    removeme = s<0.5e-2;
    TT(removeme) = [];
    p(removeme) = [];
    XX(removeme) = [];
    s(removeme) = [];
    
    removeme = s>0.5-0.5e-2;
    TT(removeme) = [];
    p(removeme) = [];
    XX(removeme) = [];
    s(removeme) = [];

    % Adding Noise
    d = s + normrnd(0,noise,size(s));
    

end