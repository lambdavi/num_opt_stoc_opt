function [xk,fk, gradfk_norm, k, xseq, btseq, deltaxk_norm] = c_steepest_desc_backtrack( ...
    x0, f, gradf, kmax, tolgrad, c1, rho, btmax, gamma, tolx, Pi_X)

% Function handle for the armijo condition
farmijo = @(fk, alpha, gradfk, pk) ...
    fk + c1 * alpha * gradfk' * pk;

% Initializations
xseq = zeros(length(x0), kmax);
btseq = zeros(1, kmax);
xk = Pi_X(x0);
fk = f(xk);
gradfk = gradf(xk);
k = 0;
gradfk_norm = norm(gradfk);
deltaxk_norm = tolx+1;

while k < kmax && gradfk_norm >= tolgrad && deltaxk_norm >= tolx
    % Compute the descent direction
    pk = -gradf(xk);
    xk_hat = Pi_X(xk+ gamma*pk);
    pi_k = xk_hat - xk;
    % Reset the value of alpha
    alpha = 1;
    
    % Compute the candidate new xk
    xnew = xk + alpha * pi_k;
    % Compute the value of f in the candidate new xk
    fnew = f(xnew);
    
    bt = 0;
    % Backtracking strategy: 
    % 2nd condition is the Armijo condition not satisfied
    while bt < btmax && fnew > farmijo(fk, alpha, gradfk, pi_k)
        % Reduce the value of alpha
        alpha = rho * alpha;
        % Update xnew and fnew w.r.t. the reduced alpha
        xnew = xk + alpha * pi_k;
        fnew = f(xnew);
        
        % Increase the counter by one
        bt = bt + 1;
        
    end

    deltaxk_norm = norm(xnew - xk);
    % Update xk, fk, gradfk_norm
    xk = xnew;
    fk = fnew;
    gradfk = gradf(xk);
    gradfk_norm = norm(gradfk);

    % Increase the step by one
    k = k + 1;
    
    % Store current xk in xseq
    xseq(:, k) = xk;
    % Store bt iterations in btseq
    btseq(k) = bt;
end

% "Cut" xseq and btseq to the correct size
xseq = xseq(:, 1:k);
btseq = btseq(1:k);

end