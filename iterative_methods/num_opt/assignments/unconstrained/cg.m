function p = cg(x0, max_iter, tol)
    r_0 = b - A*x0;
    d_k = r_0;
    x_k = x0;
    r_k=r_0;
    for i = 1:max_iter
        er = norm(r_k)/norm(r_0);
        if er <= tol
            break
        end
        z_k = A*d_k;
        alpha_k = (r_k'*d_k)/(d_k'*z_k);
        x_k = x_k + alpha_k*d_k;
        r_k1=r_k - alpha_k*z_k;
        beta_k = (r_k1'*r_k1)/(r_k'*r_k);
        d_k=r_k1+beta_k*d_k;
        r_k=r_k1;
    end
    return 
end