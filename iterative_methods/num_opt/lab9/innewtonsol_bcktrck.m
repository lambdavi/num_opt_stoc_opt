function [xk, merit_fk, merit_norm, k, xseq, btseq] = innewtonsol_bcktrck(x0, F, JF, ...
    FDJ, h, forcing_terms, gmres_maxiters_inner, gmres_maxiters_outer, kmax, tolgrad, c1, rho, btmax)
% NEWTON METHOD WITH BACKTRACKING

switch FDJ
    case 'fw'
        JF = @(x) fin_diff_J(F, x, h, 'fw');
        grad_merit = @(x) JF(x)'*F(x);

    case 'c'   
        JF = @(x) fin_diff_J(F, x, h, 'c');
        grad_merit = @(x) JF(x)'*F(x);

    case 'MF'
        JF_vector = @(x, p) (F(x + h * p) - F(x)) / h;

    otherwise
        grad_merit = @(x) JF(x)'*F(x);

end

% Defining merit functions and his gradient with function handles
merit_f= @(x)(0.5)*sum(F(x).^2);


% Function handle for the armijo condition
% Use merit function
farmijo = @(merit_fk, alpha, grad_delta) ...
    merit_fk + c1 * alpha * grad_delta;


xseq = zeros(length(x0), kmax);
btseq = zeros(1, kmax);

xk = x0;
merit_fk = merit_f(xk);
%gradmeritk = grad_merit(xk);
k = 0;
merit_norm = norm(merit_fk);
Fk=F(xk);
while k < kmax && merit_norm >= tolgrad
    tol = forcing_terms(gradfk);
    if isequal("MF", FDJ)
        J_pk = @(p) JF_vector(xk, p);
        deltaxk = gmres(J_pk, -Fk, gmres_maxiters_inner, tol, gmres_maxiters_outer);
        gradf_delta = Fk'*JF_vector(xk, deltaxk);
    else
        deltaxk = gmres(JF(xk),-Fk, gmres_maxiters_inner, tol, gmres_maxiters_outer);
        gradf_delta = grad_merit(xk)'*JF_vector(xk, deltaxk);
    end
    alpha=1;
    xnew = xk + alpha*deltaxk;
    merit_fk_new = merit_f(xnew);
    bt = 0;
    farmijo_val = farmijo(merit_fk, alpha, gradf_delta);
    while bt < btmax && merit_fk_new > farmijo_val
        % Reduce the value of alpha
        alpha = rho * alpha;
        % Update xnew and fnew w.r.t. the reduced alpha
        xnew = xk + alpha * deltaxk;
        merit_fk_new = merit_f(xnew);
        % Increase the counter by one
        bt = bt + 1;
        
    end
    xk=xnew;
    merit_fk=merit_fk_new;
    Fk=F(xk);
    merit_norm = norm(F(xk));
    k = k + 1;
    
    % Store current xk in xseq
    xseq(:, k) = xk;
    btseq(k) = bt;

end
xseq = xseq(:, 1:k);
btseq = btseq(1:k);
end