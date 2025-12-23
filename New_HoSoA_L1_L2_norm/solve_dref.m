function d0 = solve_dref(dG_dy,x_vals,d_chem,T_ref)

    options = optimoptions('fsolve','FunValCheck', 'off','MaxFunctionEvaluations',1000000,...
        'MaxIterations',1000000,'Display','off');%,'FinDiffRelStep',1e-12);
    

    d0 = zeros(length(x_vals),1);

%     for i=1:length(x_vals)
% %         DREF = matlabFunction(simplify(subs(dG_dy,x,x_vals(i))));
% %         coords = d_chem + subs(dG_dy,x,x_vals(i)) == 0;
%         coords = @(z) d_chem + dG_dy(T_ref,x_vals(i),z);
%         dref = fsolve(coords,X0(i),options);
% %         dref = solve(coords,z);      
%   
%         if double(dref)<0 || ~isreal(double(dref))
%             d0(i) = 0;
%         else
%             d0(i) = dref;
%         end
%     end    
    
    for i=1:length(x_vals)
%         xo = X0(i);
        %dG_dd = subs(sym(dG_dd_og),x,x_vals(i));
        %dG_dd = matlabFunction(expand(dG_dd));
        coords = @(z) d_chem + dG_dy(T_ref,x_vals(i),z);
        %dref = fsolve(coords,xo,optimset('FunValCheck', 'off', 'Display', 'off'));
        if coords(1e-12) < 0 && coords(0.5 - 1e-12) > 0
            dref = fzero(coords,[1e-12 0.5-1e-12],optimset('Display', 'off'));
            if dref < 1e-3
                d0(i) = 0.001;
            else
                d0(i) = dref;
            end
        else
            d0(i) = 0.001;
        end 
    end
   
end

