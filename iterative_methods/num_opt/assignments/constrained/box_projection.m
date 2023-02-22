function [xhat] = box_projection(x, mins, maxs)
    xhat = min(maxs, max(x, mins));
end

