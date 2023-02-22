function [xk, fk, gradfk_norm, k, xseq, btseq, gmres_it, ratio] = newton_bcktrck(x0, f, gradf, Hessf, ...
    kmax, tolgrad, c1, rho, btmax,mem)
% NEWTON METHOD WITH BACKTRACKING


% Function handle for the armijo condition
farmijo = @(fk, alpha, gradfk, pk) ...
    fk + c1 * alpha * gradfk' * pk;

if mem == 1
        xseq = zeros(length(x0), kmax);
else
    xseq=0;
end
btseq = zeros(1, kmax);
gmres_it = zeros(1, kmax);
xk = x0;
fk = zeros(1, kmax+1);
fk(1) = f(xk);
k = 0;
gradfk = gradf(xk);
gradfk_norm = zeros(kmax+1, 1); 
ratio = zeros(kmax-4,1);
gradfk_norm(1) = norm(gradfk);
while k < kmax && gradfk_norm(k+1) >= tolgrad

    H=Hessf(xk);
    if class(H) == "sparse"
        L = ichol(H,struct('diagcomp',0.01)); % create a diagonal preconditioner
        [pk,flag,~,iter] = pcg(H, -gradfk, 1e-6, 50,L,L');
    else
        [pk,flag,~,iter] = pcg(H, -gradfk, 1e-6, 50);
    end
    
    alpha=1;
    xnew = xk + alpha*pk;
    fnew = f(xnew);
    bt = 0;
    while bt < btmax && fnew > farmijo(fk(k+1), alpha, gradfk, pk)
        % Reduce the value of alpha
        alpha = rho * alpha;
        % Update xnew and fnew w.r.t. the reduced alpha
        xnew = xk + alpha * pk;
        fnew = f(xnew);
        % Increase the counter by one
        bt = bt + 1;
        
    end
    k = k + 1;

    xk=xnew;
    fk(k+1)=fnew;
    if k >= 4
        ratio(k-3) = (fk(k-2) - fk(k-1)) / (fk(k-1) - fk(k));
    end
    gradfk = gradf(xk);
    gradfk_norm(k+1) = norm(gradfk);
    
    if mem == 1
            % Store current xk in xseq
            xseq(:, k) = xk;
    end
        % Store bt iterations in btseq
        btseq(k) = bt;
        gmres_it(k) = iter;


if mem == 1
        xseq = xseq(:, 1:k);
end
btseq = btseq(1:k);
gradfk_norm = gradfk_norm(1:k+1);
gmres_it = gmres_it(1:k);
fk = fk(1:k+1);
ratio = ratio(1:k-3);
end