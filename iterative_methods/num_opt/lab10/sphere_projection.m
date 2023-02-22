function [xhat] = sphere_projection(x,c,r)
    norm_xc = norm(x-c);
    if norm_xc <= r
        xhat = x;
    else
        xhat = c + r*(x-c)/norm_xc;
    end
   
end

