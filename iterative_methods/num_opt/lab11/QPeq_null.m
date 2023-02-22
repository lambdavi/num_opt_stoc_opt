function [xstar,fxstar, lambda_star, v_star, KKT_gradL_norm, KKT_eq_norm] = QPeq_null( ...
    Q,A,x2)
    [r,c] = size(A);
    A1=A(:,1:r);
    A2=A(:,r+1:end);
    
    Z = [-A1\A2; eye(n-m)];

    %init xhat
    xhat = []
    % compute v_star

    % compite x_star given vstar and xhat

    %compute fxstar ecc
end

