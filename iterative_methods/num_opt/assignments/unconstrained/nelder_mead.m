function [x_min, f_min] = nelder_mead(func, x0, tol, max_iter)
    n = length(x0);
    simplex = zeros(n + 1, n);
    simplex(1, :) = x0;
    for i = 2 : n + 1
        eye_matrix = eye(n);
        simplex(i, :) = x0 + eye_matrix(i - 1, :);

    end
    for i = 1 : max_iter
        [~, indices] = sort(arrayfun(func, simplex));
        simplex = simplex(indices, :);
        x_bar = mean(simplex(1 : n, :), 1);
        x_r = 2 * x_bar - simplex(n + 1, :);
        f_r = func(x_r);
        if f_r < func(simplex(1, :))
            x_e = 3 * x_bar - 2 * simplex(n + 1, :);
            f_e = func(x_e);
            if f_e < f_r
                simplex(n + 1, :) = x_e;
            else
                simplex(n + 1, :) = x_r;
            end
        else
            if f_r < func(simplex(n, :))
                simplex(n + 1, :) = x_r;
            else
                x_c = 0.5 * x_bar + 0.5 * simplex(n + 1, :);
                f_c = func(x_c);
                if f_c < func(simplex(n + 1, :))
                    simplex(n + 1, :) = x_c;
                else
                    simplex(2 : n + 1, :) = 0.5 * (simplex(2 : n + 1, :) + simplex(ones(n, 1), :));
                end
            end
        end
        if abs(func(simplex(n + 1, :)) - func(simplex(1, :))) < tol
            break;
        end
    end
    x_min = simplex(1, :);
    f_min = func(x_min);
end
