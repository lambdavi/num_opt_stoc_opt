clear
clc

n=50;
m=10;
A=rand(m,n);
b=rand(m,1);
c=rand(n,1);
Q=diag(diag(rand(n)));


[xstar, fstar, lambda_star, KKT_gradL_norm, KKT_eq_norm] = QPeq_Schur(Q, inv(Q), c, A, b);

[xstar_ai, fstar_ai, lambda_star_ai, KKT_gradL_norm_ai, KKT_eq_norm_ai] = QPeq_Schur_autoinv(Q, c, A, b);