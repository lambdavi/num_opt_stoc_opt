n = 1000;
alpha = 2;
tol = 1e-4;
a = [-ones(n,1) alpha*ones(n,1) -ones(n,1)];
A = spdiags(a,-1:1,n,n);
b = A*ones(n,1);
x_0 = zeros(n,1);

% PCG Implementation
[x, flag, relres, iter] = pcg(A,b,tol);
disp(iter)

% "From zero" implementation
r_0 = b - A*x_0;
d_k = r_0;
x_k = x_0;
r_k=r_0;
for i = 1:10000
    er = norm(r_k)/norm(r_0);
    %disp(er)
    if er <=tol
        break
    end
    z_k = A*d_k;
    alpha_k = (r_k'*d_k)/(d_k'*z_k);
    x_k = x_k + alpha_k*d_k;
    r_k1=r_k - alpha_k*z_k;
    beta_k = (r_k1'*r_k1)/(r_k'*r_k);
    d_k=r_k1+beta_k*d_k;
    r_k=r_k1;
    fk = f(x_k);
end
disp(i)


