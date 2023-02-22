function [JFx] = fin_diff_J(F, x, h, type)
%
% function [gradf] = findiff_grad(f, x, h, type)
%
% Function that approximate the gradient of f in x (column vector) with the
% finite difference (forward/centered) method.
%
% INPUTS:
% f = function handle that describes a function R^n->R;
% x = n-dimensional column vector;
% h = the h used for the finite difference computation of gradf
% type = 'fw' or 'c' for choosing the forward/centered finite difference
% computation of the gradient.
%
% OUTPUTS:
% gradfx = column vector (same size of x) corresponding to the approximation
% of the gradient of f in x.

Fx = F(x);
JFx = zeros(length(Fx), length(x));

switch type
    case 'fw'
        for j=1:length(Fx)
            xh = x;
            xh(j) = xh(j) + h;
            JFx(:,j) = (F(xh) - Fx)/ h;
        end
    case 'c'
        for j=1:length(Fx)
            xh_plus = x;
            xh_minus = x;
            xh_plus(j) = xh_plus(j) + h;
            xh_minus(j) = xh_minus(j) - h;
            JFx(:, j)=(F(xh_plus) - F(xh_minus))/(2*h);
        end
    otherwise % repeat the 'fw' case
        disp("Not specified the method fw/c")
 end
end

