function [xk, fk, gradfk_norm, k, xseq, btseq] = mod_newton_bcktrck(x0, f, gradf, Hessf, ...
    kmax, tolgrad, c1, rho, btmax)
% MODIFIED NEWTON METHOD WITH BACKTRACKING


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
    H = Hessf(xk);
    pk = -H\gradfk;
    while (pk' * gradfk) > 0
        beta = norm(H,"fro");
        % check the diagonal for minimum
        min = 9999;
        for i=length(H)
            for j = length(H)
                if H(i,j) < min
                    min = H(i,j);
                end
            end
        end
        if min > 0
            tau_k=0;
        else
            tau_k = beta/2;
        end

        B = H + tau_k*ones(size(H));
        [R,flag] = chol(B);
        if flag == 0
            break;
        end
         tau_k = max(2*tau_k, beta/2);
    end
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
end
