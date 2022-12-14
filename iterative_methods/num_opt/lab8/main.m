load("test_nonlinsys.mat")

x0 = [-456, 1000]';
c1=1e-4;
rho=0.8;
btmax=50;
tolgrad = 1e-12;
[xk, fk, gradfk_norm, k, xseq, btseq] = newtonsol_bcktrck(x0, F, JF, ...
    kmax, tolgrad, c1, rho, btmax);

