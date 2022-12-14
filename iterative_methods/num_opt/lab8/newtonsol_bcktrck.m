function [xk, merit_fk, merit_norm, k, xseq, btseq] = newtonsol_bcktrck(x0, F, JF, ...
    kmax, tolgrad, c1, rho, btmax)
% NEWTON METHOD WITH BACKTRACKING


% Defining merit functions and his gradient with function handles
merit_f= @(x)(0.5)*sum(F(x).^2);
grad_merit = @(x) JF(x)'*F(x);

% Function handle for the armijo condition
% Use merit function
farmijo = @(merit_fk, alpha, grad_meritk, deltaxk) ...
    merit_fk + c1 * alpha * grad_meritk' * deltaxk;


xseq = zeros(length(x0), kmax);
btseq = zeros(1, kmax);

xk = x0;
merit_fk = merit_f(xk);
gradmeritk = grad_merit(xk);
k = 0;
merit_norm = norm(merit_fk);
while k < kmax && merit_norm >= tolgrad
    deltaxk = -JF(xk)\F(xk);
    alpha=1;
    xnew = xk + alpha*deltaxk;
    merit_fk_new = merit_f(xnew);
    bt = 0;
    farmijo_val = farmijo(merit_fk, alpha, gradmeritk, deltaxk);
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
    gradmeritk = grad_merit(xk);
    merit_norm = norm(F(xk));
    k = k + 1;
    
    % Store current xk in xseq
    xseq(:, k) = xk;
    btseq(k) = bt;

end
xseq = xseq(:, 1:k);
btseq = btseq(1:k);
end