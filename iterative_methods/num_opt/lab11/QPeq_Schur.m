function [xstar, fstar, lambda_star, KKT_gradL_norm, KKT_eq_norm] = QPeq_Schur( ...
    Q, Qinv, c, A, b)
    S = A*Qinv*A';
    beta = -b -A*Qinv*c;
    lambda_star = S\beta;
    xstar = Qinv*(-c -A'*lambda_star);
    fstar = 0.5*xstar'*Q*xstar + c'*xstar;
    KKT_gradL_norm = norm(Q*xstar + c+A'*lambda_star);
    KKT_eq_norm = norm(A*xstar-b);
end

