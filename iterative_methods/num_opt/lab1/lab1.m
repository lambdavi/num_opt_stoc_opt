n = 1000;
alpha = 4;
tol = 1e-4;
a = [-ones(n,1) alpha*ones(n,1) -ones(n,1)];
A = spdiags(a,-1:1,n,n);
b = A*ones(n,1);
x_0 = zeros(n,1);
r_0 = b-A*x_0;
r_k = r_0;
x_k = x_0;
er = 100;
for i = 1:10000
    er = norm(r_k)/norm(r_0);
    disp(er)
    if er <=tol
        break
    end
    z_k = A*r_k;
    alpha_k = (r_k'*r_k)/(r_k'*z_k);
    x_k = x_k + alpha_k * r_k;
    r_k = r_k - alpha_k * z_k;
end
disp('Number of steps: ')
disp(i)
