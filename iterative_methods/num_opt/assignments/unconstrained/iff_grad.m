function grad = iff_grad(x, h, diff_type)
    %h=h*norm(x);
    [n, ~] = size(x);
    seq = (1:n)';
    switch diff_type
        case 'fw'
            grad = seq .* ((x+h).^2 - (x).^2) ./ h;
        case 'c'
            grad = seq .* ((x+h).^2 - (x-h).^2) ./ (2.*h);
        otherwise
            error('Requested difference not available')
    end
end