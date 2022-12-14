function [Hessfx] = fin_diff_Hess(f, x, h)
%
n = length(x);
Hessfx = zeros(n);
h = sqrt(h);
% Diagonals
for c=1:n
    xk_p = x;
    xk_m = x;
    xk_p(c) = x(c)+h;
    xk_m(c) = x(c)-h;
    Hessfx(c,c) = (f(xk_p) - 2*f(x) + f(xk_m))/(h^2);
end

for c=1:n
    r=1;
    while r<c
        xk_r = x;
        xk_c = x;
        xk_both = x;

        xk_both(c) = x(c)+h;
        xk_both(r) = x(r)+h;

        xk_r(r) = x(r)+h;
        xk_c(c) = x(c)+h;

        Hessfx(c,r) = (f(xk_both) - f(xk_r) - f(xk_c) + f(x))/(h^2);
        Hessfx(r,c) = Hessfx(c,r);
        r=r+1;
    end
end