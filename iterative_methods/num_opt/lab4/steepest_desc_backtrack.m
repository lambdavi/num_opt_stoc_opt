function [xk,fk, gradfk_norm, k, xseq, btseq] = steepest_desc_backtrack(x0, f, gradf, alpha0, kmax, tolgrad, c1, rho, btmax)
% STEEPEST DESCENT METHOD WITH BACKTRACKING 
% that implements the steepest descent optimization method with the backtracking strategy
% INPUTS:
% x0: a column vector of n elements representing the starting point for the optimization method;
% f: a function handle variable that, for each column vector x ∈ R n, returns the value f(x), where f : R^n → R is the loss function that have to be minimized;
% gradf: a function handle variable that, for each column vector x ∈ R^n, returns the value ∇f(x) as a column vector, where ∇f : R^n → R^n is the gradient of f;
% alpha0: a real scalar value characterizing the step length of the optimization method;
% kmax: an integer scalar value characterizing the maximum number of iterations of the method;
% tolgrad: a real scalar value characterizing the tolerance with respect to the norm of the gradient in order to stop the method.
% c1: the factor c1 for the Armijo condition that must be a scalar in (0, 1);
% rho: fixed factor, less than 1, used to reduce α;
% btmax: maximum number of steps allowed to update α during the backtracking strategy
%

% Initializations
xseq = zeros(length(x0), kmax);
btseq = zeros(1,btmax);
xk = x0;
alphak=alpha0;
k = 0;
gradfk_norm = norm(gradf(xk));

while k < kmax && gradfk_norm >= tolgrad
    % Compute the descent direction
    pk = -gradf(xk);
    j=0;
    alphak=alpha0;
    % Backtracking logic
    while j <= btmax
        % Compute the new value for xk
        xk_1 = xk + alphak * pk;
        % Armijo condition
        if f(xk_1) <= f(xk) + c1*alphak*gradf(xk)'*pk
            disp(alphak)
            xk=xk_1;
            break;
        else
            % Update of alphak using rho
            alphak=rho*alphak;
            j=j+1;
        end
    end
    btseq(k+1) = j; 

    % Compute the gradient of f in xk
    gradfk_norm = norm(gradf(xk));
    
    % Increase the step by one
    k = k + 1;
    
    % Store current xk in xseq
    xseq(:, k) = xk;
end

% Compute f(xk)
fk = f(xk);

% "Cut" xseq to the correct size
xseq = xseq(:, 1:k);
end