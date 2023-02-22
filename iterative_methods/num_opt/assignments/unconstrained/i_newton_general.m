function [xk, fk, gradfk_norm, k, xseq, btseq] = ...
    i_newton_general(x0, f, gradf, Hessf, kmax, ...
    tolgrad, c1, rho, btmax, FDgrad, FDHess, h, pcg_maxit, forcing_terms)
%
%
% [xk, fk, gradfk_norm, k, xseq] = ...
%     newton_general(x0, f, gradf, Hessf, ...
%     kmax, tolgrad, c1, rho, btmax, FDgrad, FDhess, h, pcg_maxit)
%
% Function that performs the newton optimization method, 
% implementing the backtracking strategy and, optionally, finite
% differences approximations for the gradient and/or the Hessian.
%
% INPUTS:
% x0 = n-dimensional column vector;
% f = function handle that describes a function R^n->R;
% gradf = function handle that describes the gradient of f (not necessarily
% used);
% Hessf = function handle that describes the Hessian of f (not necessarily 
% used);
% kmax = maximum number of iterations permitted;
% tolgrad = value used as stopping criterion w.r.t. the norm of the
% gradient;
% c1 = ﻿the factor of the Armijo condition that must be a scalar in (0,1);
% rho = ﻿fixed factor, lesser than 1, used for reducing alpha0;
% btmax = ﻿maximum number of steps for updating alpha during the 
% backtracking strategy;
% FDgrad = 'fw' (FD Forward approx. for gradf), 'c' (FD Centered approx. 
% for gradf), any other string (usage of input Hessf)
% FDHess = 'fw' (FD approx. for Hessf), 'Jfw' (Jacobian FD Forward
% approx. of Hessf), 'Jc' (Jacobian FD Centered approx. of Hessf), 'MF'
% (Matrix Free implementation for solving Hessf(xk)pk=-gradf(xk)), any 
% other string (usage of input Hessf);
% h = approximation step for FD (if used);
% pcg_maxit = maximum number of iterations for the pcg solver.
%
% OUTPUTS:
% xk = the last x computed by the function;
% fk = the value f(xk);
% gradfk_norm = value of the norm of gradf(xk)
% k = index of the last iteration performed
% xseq = n-by-k matrix where the columns are the xk computed during the 
% iterations
% btseq = 1-by-k vector where elements are the number of backtracking
% iterations at each optimization step.
%
switch FDgrad
    case 'fw'
        % OVERWRITE gradf WITH A F. HANDLE THAT USES findiff_grad
        % (with option 'fw')
        gradf = @(x) fin_diff_grad(f, x, h, 'fw');
        
    case 'c'
        % OVERWRITE gradf WITH A F. HANDLE THAT USES findiff_grad
        % (with option 'c')        
        gradf = @(x) fin_diff_grad(f, x, h, 'c');
        
    otherwise
        % WE USE THE INPUT FUNCTION HANDLE gradf...
        %
        % THEN WE DO NOT NEED TO WRITE ANYTHING!
        % ACTUALLY WE COULD DELETE THE OTHERWISE BLOCK
        
end

% (OPTIONAL): IN CASE OF APPROXIMATED GRADIENT, IT IS BETTER TO NOT
% APPROXIMATE Hessf WITH THE JACOBIAN!
if isequal(FDgrad, 'fw') || isequal(FDgrad, 'c')
    switch FDHess
        case 'fw'
            % OVERWRITE Hessf WITH A F. HANDLE THAT USES findiff_Hess
            Hessf = @(x) fin_diff_Hess(f, x, sqrt(h));
        case 'MF'
            % DEFINE a f. handle for the product of Hessf * p USING THE
            % GRADIENT
            Hessf_pk = @(x, p) (gradf(x + h * p) - gradf(x)) / h;
        otherwise
            % WE USE THE INPUT FUNCTION HANDLE Hessf
            %
            % THEN WE DO NOT NEED TO WRITE ANYTHING!
            % ACTUALLY WE COULD DELETE THE OTHERWISE BLOCK
    end
else
    switch FDHess
        case 'fw'
            % OVERWRITE Hessf WITH A F. HANDLE THAT USES findiff_Hess
            Hessf = @(x) fin_diff_Hess(f, x, sqrt(h));
        case 'Jfw'
            % OVERWRITE Hessf WITH A F. HANDLE THAT USES findiff_J
            % (with option 'fw')
            Hessf = @(x) fin_diff_J(gradf, x, h, 'fw');

        case 'Jc'
            % OVERWRITE Hessf WITH A F. HANDLE THAT USES findiff_J
            % (with option 'c')
            Hessf = @(x) findiff_J(gradf, x, h, 'c');

        case 'MF'
            % DEFINE a f. handle for the product of Hessf * p USING THE
            % GRADIENT
            Hessf_pk = @(x, p) (gradf(x + h * p) - gradf(x)) / h;

        otherwise
            % WE USE THE INPUT FUNCTION HANDLE Hessf
            %
            % THEN WE DO NOT NEED TO WRITE ANYTHING!
            % ACTUALLY WE COULD DELETE THE OTHERWISE BLOCK
    end
end

% Function handle for the armijo condition
farmijo = @(fk, alpha, gradfk, pk) ...
    fk + c1 * alpha * gradfk' * pk;

% Initializations
xseq = zeros(length(x0), kmax);
btseq = zeros(1, kmax);

xk = x0;
fk = zeros(1, kmax+1);
fk(1) = f(xk);

k = 0;
gradfk = gradf(xk);
gradfk_norm = zeros(kmax+1, 1); 
gradfk_norm(1) = norm(gradfk);


while k < kmax && gradfk_norm(k+1) >= tolgrad
    % Compute the descent direction as solution of
    % Hessf(xk) p = - graf(xk)
    tol_pcg = forcing_terms(gradfk);
    switch FDHess
        case 'MF'
            % ITERATIVE METHOD (DIRECTED METHODS DO NOT WORK WITH MF 
            % IMPLEMENTATION)
            % OBSERVATION: We use pcg with a f. handle as first argument
            % that describes the linear product Hessf(xk) p.
            %
            % Then, we define a new f. handle (to read better the code)
            % that exploits the f. handle Hessf_pk to define the product
            % between Hessfk and p (i.e., xk fixed, p variable):
            Hessfk_pk = @(p) Hessf_pk(xk, p);
            [pk,FLAG,RELRES,ITER] = pcg(Hessfk_pk, -gradfk, tol_pcg, pcg_maxit);
            
            % ALTERNATIVE (ONE LINE OF CODE):
            % pk = pcg(@(p)Hessf_pk(xk, p), -gradfk, tolgrad, pcg_maxit);

        otherwise
            % DIRECT METHOD (uncomment to use this method)
            % pk = -Hessf(xk)\gradf(xk)            
            
            % ITERATIVE METHOD (uncomment to use this method)
            % OBERVATION: simple usage of pcg with matrix of the linear
            % system as first input of the pcg function.
            H=Hessf(xk);
            L = ichol(H,struct('diagcomp',0.01)); % create a diagonal preconditioner
            [pk,FLAG,RELRES,ITER] = pcg(H, -gradfk, tol_pcg, pcg_maxit,L,L');
    end
    
    
    % Reset the value of alpha
    alpha = 1;
    
    % Compute the candidate new xk
    xnew = xk + alpha * pk;
    % Compute the value of f in the candidate new xk
    fnew = f(xnew);
    
    bt = 0;
    % Backtracking strategy: 
    % 2nd condition is the Armijo condition not satisfied
    while bt < btmax && fnew > farmijo(fk(k+1), alpha, xk, pk)
        % Reduce the value of alpha
        alpha = rho * alpha;
        % Update xnew and fnew w.r.t. the reduced alpha
        xnew = xk + alpha * pk;
        fnew = f(xnew);
        
        % Increase the counter by one
        bt = bt + 1;
        
    end
    k = k + 1;

    % Update xk, fk, gradfk_norm
    xk = xnew;
    fk(k+1)=fnew;
    gradfk = gradf(xk);
    gradfk_norm(k+1) = norm(gradfk);    
    % Increase the step by one
    
    % Store current xk in xseq
    xseq(:, k) = xk;
    % Store bt iterations in btseq
    btseq(k) = bt;
end

% "Cut" xseq and btseq to the correct size
xseq = xseq(:, 1:k);
btseq = btseq(1:k);
gradfk_norm = gradfk_norm(1:k+1);
fk = fk(1:k+1);



end