function [gradfx] = fin_diff_grad(f, x, h, type)
%FIND_IFF_GRAD 
n = length(x);
gradfx = zeros(n,1);
switch type
    case 'fw'
        for i=1:n
            xk = x;
            xk(i) = xk(i)+h;
            gradfx(i) = (f(xk)- f(x))/h;
        end
    case 'c'
        for i=1:n
            xk_p = x;
            xk_m = x;
            xk_p(i) = x(i)+h;
            xk_m(i) = x(i)-h;
            gradfx(i) = (f(xk_p)- f(xk_m))/(2*h);
        end
    otherwise
        disp("Error in type field");
end

end

