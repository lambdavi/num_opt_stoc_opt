function [xk, fk, gradfk_norm, k, xseq, btseq] = newton_general(x0, f, gradf, Hessf, ...
    kmax, tolgrad, c1, rho, btmax, FDgrad, FDHess, h)
% NEWTON METHOD WITH BACKTRACKING

switch FDgrad 
    case 'fw'
        % OVERWRITE THE FUNCTION HANDLE GRADF 
        gradf = @(x) fin_diff_grad(f,x,h,FDgrad);
    case 'c'
        gradf = @(x) fin_diff_grad(f,x,h,FDgrad);
    otherwise
        disp("Using standard gradf")
end

switch FDHess
    case 'fw'
        % OVERWRITE THE FUNCTION HANDLE GRADF 
        Hessf = @(x) fin_diff_Hess(f,x,sqrt(h));
    otherwise
        disp("Using standard Hess")
end
% Function handle for the armijo condition
farmijo = @(fk, alpha, gradfk, pk) ...
    fk + c1 * alpha * gradfk' * pk;


xseq = zeros(length(x0), kmax);
btseq = zeros(1, kmax);

xk = x0;
fk = f(xk);
gradfk = gradf(xk);
k = 0;
gradfk_norm = norm(gradfk);
while k < kmax && gradfk_norm >= tolgrad
    pk = -Hessf(xk)\gradfk;
    alpha=1;
    xnew = xk + alpha*pk;
    fnew = f(xnew);
    bt = 0;
    while bt < btmax && fnew > farmijo(fk, alpha, gradfk, pk)
        % Reduce the value of alpha
        alpha = rho * alpha;
        % Update xnew and fnew w.r.t. the reduced alpha
        xnew = xk + alpha * pk;
        fnew = f(xnew);
        % Increase the counter by one
        bt = bt + 1;
        
    end
    xk=xnew;
    fk=fnew;
    gradfk = gradf(xk);
    gradfk_norm = norm(gradfk);
    k = k + 1;
    
    % Store current xk in xseq
    xseq(:, k) = xk;
    btseq(k) = bt;

end
xseq = xseq(:, 1:k);
btseq = btseq(1:k);
end