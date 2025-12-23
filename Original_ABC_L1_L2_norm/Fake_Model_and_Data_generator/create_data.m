function [TT, p, XX, YY, d, s] = create_data(G,N,noise,doop)

    syms x y z T muh real
    
    % chempot = -diff(G,z);
    % chempot = matlabFunction(chempot);

    chempot = -diff(G,z) - muh;
    % eqn = solve(chempot,z);
    solMu = matlabFunction(chempot);

    T1 = linspace(800,1800,N);
    % delta = linspace(0.01,0.3,N);
    X = linspace(1e-6,1-1e-6,N);
    Y = linspace(1e-6,1-1e-6,N);
    % P = logspace(-4,-0.6778,N);
    P = logspace(-10,-.6777777,N);
    
    % [TT, DD, XX] = meshgrid(T1,delta,X);
    [TT, PP, XX, YY] = ndgrid(T1,P,X,Y);

    TT = reshape(TT,[1 N^4]);
    % DD = reshape(DD,[1 N^3]);
    XX = reshape(XX,[1 N^4]);
    YY = reshape(YY,[1 N^4]);
    PP = reshape(PP,[1 N^4]);

    % make duplicates
    if doop > 1
        pp = PP;
        xx = XX;
        tt = TT;
        yy = YY;
        for i = 1:doop-1
        PP = [PP,pp];
        XX = [XX,xx];
        TT = [TT,tt];
        YY = [YY,yy];
        end
    end
    
    omuh = solve_for_muO(TT,PP);
    s = 0*omuh;
    for i = 1:length(TT)
        % eqn = subs(chempot,[T x muh],[TT(i) XX(i) omuh(i)]);
        if solMu(TT(i),omuh(i),XX(i),YY(i),1e-10) > 0 && solMu(TT(i),omuh(i),XX(i),YY(i),0.5-1e-10) < 0
            eqn = @(z) solMu(TT(i),omuh(i),XX(i),YY(i),z);
            s(i) = double(fzero(eqn,[1e-10 0.5-1e-10]));
        else
            s(i) = 0;
        end
    end
    
    p = PP;
    removeme = s<0.5e-2;
    TT(removeme) = [];
    p(removeme) = [];
    XX(removeme) = [];
    YY(removeme) = [];
    s(removeme) = [];
    
    removeme = s>0.5-0.5e-2;
    TT(removeme) = [];
    p(removeme) = [];
    XX(removeme) = [];
    YY(removeme) = [];
    s(removeme) = [];
    
    % if you want to solve for P
    % ochem = chempot(TT,XX,DD);
    % p = solve_for_p(TT,ochem);

    % Adding Noise
    d = s + normrnd(0,noise,size(s));
    
%     removeme = d<0.5e-2;
%     d(removeme) = [];
%     TT(removeme) = [];
%     p(removeme) = [];
%     XX(removeme) = [];
%     s(removeme) = [];
%     
%     removeme = d>0.5-0.5e-2;
%     d(removeme) = [];
%     TT(removeme) = [];
%     p(removeme) = [];
%     XX(removeme) = [];
%     s(removeme) = [];

end