% Function for computing the value of F(x) at point X where F(x) is the one
% reported in problem 79 in [1]
% INPUTS:
% X: n by 1 vector representing a point in the n-dimensional space
% OUTPUTS:
% Fx: the value of the function at X

x = [1;1;1;1];

Fx = @(x) 0.5*sum(small_fx_79(x).^2);
gradf = @(x) problem_79_grad(x);
Hessf = @(x) problem_79_hess(x);

function small_fx_array = small_fx_79(x)
    n = length(x);
    small_fx_array = zeros(n,1);
    for k = 1:n
        small_fx_array(k) = (3*x(k) - 0.1*x(k).^2 + 1);
        if(k==1)
            small_fx_array(k) = small_fx_array(k) - 2*x(k+1);
        elseif (k==n)
            small_fx_array(k) = small_fx_array(k) - x(k-1);
        else
            small_fx_array(k) = small_fx_array(k) - x(k-1) - 2*x(k+1);
        end
    end
end

function gradFx = problem_79_grad(x)
% Function for computing the gradient vector at a specified point of the
% function reported in problem 76 in [1]
% INPUTS:
% X: n by 1 vector representing a point in the n-dimensional space
% OUTPUTS:
% gradFx: n by 1 vector which represent the gradient at X
    small_fx_array = small_fx_79(x);
    n = length(x);
    gradFx = zeros(n,1);
    gradFx(1) = (3 - 0.2*x(1))*small_fx_array(1) - small_fx_array(2);
    gradFx(n) = (3 - 0.2*x(n))*small_fx_array(n) - 2*small_fx_array(n-1);
    for i = 2:(n-1)
        gradFx(i) = -2*small_fx_array(i-1) + (3 - 0.2*x(i))*small_fx_array(i) - small_fx_array(i+1);
    end 
end

function HessFx = problem_79_hess(x)
% Function for computing the gradient vector at a specified point of the
% function reported in problem 76 in [1]
% INPUTS:
% X: n by 1 vector representing a point in the n-dimensional space
% OUTPUTS:
% gradFx: n by 1 vector which represent the gradient at X
    small_fx_array = small_fx_79(x);
    n = length(x);
    
    %Initialize vectors for sparse matrix:
    rows=zeros(5*n-6,1);
    cols=zeros(5*n-6,1);
    vals=zeros(5*n-6,1);
    count=3;
    % Compute the diagonal first and last element
    % HessFx(1,1) = ((3 - 0.2*x(1)).^2)*(-0.2*small_fx_array(1)) + 1;
    rows(1) = 1; cols(1) = 1; vals(1) = ((3 - 0.2*x(1)).^2)*(-0.2*small_fx_array(1)) + 1;
    % HessFx(n,n) = 4 + ((3 - 0.2*x(n)).^2)*(-0.2*small_fx_array(n));
    rows(2) = n; cols(2) = n; vals(2) = 4 + ((3 - 0.2*x(n)).^2)*(-0.2*small_fx_array(n));
    for i = 1:(n-1)
        if(i~=1)
            % Compute the diagonal
            % HessFx(i,i) = 5 + ((3 - 0.2*x(i)).^2)*(-0.2*small_fx_array(i));
            rows(count) = i; cols(count) = i; vals(count) = 5 + ((3 - 0.2*x(i)).^2)*(-0.2*small_fx_array(i));
            count = count+1;
        end
        % Compute the "second" diagonal
        % HessFx(i,i+1) = -2*(3 - 0.2*x(i))-3 - 0.2*x(i+1);
        rows(count) = i; cols(count) = i+1; vals(count) = -2*(3 - 0.2*x(i))-3 - 0.2*x(i+1);
        count = count+1;
        % HessFx(i+1,i) = HessFx(i,i+1);
        rows(count) = i+1; cols(count) = i; vals(count) = -2*(3 - 0.2*x(i))-3 - 0.2*x(i+1);
        count = count+1;

        % Compute the "third" diagonal
        if(i~=n-1)
            % HessFx(i,i+2) = 2;
            rows(count) = i; cols(count) = i+2; vals(count) = 2;
            count = count+1;
            % HessFx(i+2,i) = 2;
            rows(count) = i+2; cols(count) = i; vals(count) = 2;
            count = count+1;
        end
    end 
    HessFx = sparse(rows,cols,vals,n,n);
end
