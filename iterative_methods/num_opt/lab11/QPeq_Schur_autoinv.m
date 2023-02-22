function [xstar, fstar, lambda_star, KKT_gradL_norm, KKT_eq_norm] = QPeq_Schur_autoinv( ...
    Q, c, A, b)
    
    QAt_inv = (Q\A');
    Qc_inv = (Q\c);

    S = A * QAt_inv;
    beta = -b -A* Qc_inv;
    lambda_star = S\beta;
    xstar = -Qc_inv -QAt_inv*lambda_star;
    fstar = 0.5*xstar'*Q*xstar + c'*xstar;
    KKT_gradL_norm = norm(Q*xstar + c+A'*lambda_star);
    KKT_eq_norm = norm(A*xstar-b);
end


