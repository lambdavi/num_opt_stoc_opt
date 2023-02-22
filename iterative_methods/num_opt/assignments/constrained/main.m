clear
close all

load("test_functions_fin.mat")
c1=1e-4;
rho=0.5;
btmax=50;
gamma=1;
tolx=1e-8;
tolgrad=1e-8;
h=1e-12;
n=1e3;

mins0=-5.12*ones(n,1);
maxs0=5.12*ones(n,1);

mins1=1*ones(n,1);
maxs1=5.12*ones(n,1);

mins2=[-5.12;1*ones(n-1,1)];
maxs2=5.12*ones(n,1);

mins3=[-5.12*ones(n/2,1);1*ones(n/2,1)];
maxs3=5.12*ones(n,1);

constrm = mins1;
constrM = maxs1;

kmax=1e4;
proj = @(x) box_projection2(x,constrm, constrM);
f = @(x) weighted_sphere(x);
x0=[-2;-4*ones(n-1,1)];
gradf = @(x) grad_weighted_sphere(x);
%gradf = @(x) iff_grad_weighted_sphere(x,h,'central');

[xk, fk, gradfk_norm, deltaxk_norm, k, xseq, btseq, fseq] = ...
    constr_steepest_desc_bcktrck(x0, f, gradf, ...
    kmax, tolgrad, c1, rho, btmax, gamma, tolx, proj);

[X, Y] = meshgrid(linspace(-5.12, 5.12, 500));
Z = X.^2 + 2*Y.^2;
fig1 = figure();
% Contour plot with curve levels for each point in xseq
[C1, ~] = contour(X,Y,Z);
hold on
plot([x0(1), xseq(1,:)], [x0(2), xseq(2,:)],'--+')
rectangle('Position',[constrm(1:2)', (constrM(1:2)-constrm(1:2))' ], 'EdgeColor','r')
xlabel("X1");
ylabel("X2");

hold off



igure();
iter = 1:k; % create a vector for the number of iterations
hold on
plot(iter, fseq2(1:k))% plot the function values against the number of iterations
plot(iter, fseq(1:k))
xlabel('Number of Iterations') % label the x-axis
ylabel('Function Value') % label the y-axis
title('Function Value vs Number of Iterations') % add a title to the plot
legend('Sequence 1', 'Sequence 2')
ylim([0 50])


figure();
surf(X, Y, Z, 'EdgeColor', 'none')
plot3([x0(1) xseq(1, :)], [x0(2) xseq(2, :)], [f1(x0), f1(xseq)], 'r--*')


function y = weighted_sphere(x)
    [n, ~] = size(x);
    y = sum([1:n]' .* (x .^ 2));
end

function grad = grad_weighted_sphere(x)
    [n, ~] = size(x);
    grad = [2:2:2*n]' .* x;
end

function grad = iff_grad_weighted_sphere(x, h, diff_type)
    h=h*norm(x);
    [n, ~] = size(x);
    seq = (1:n)';
    switch diff_type
        case 'forward'
            grad = seq .* ((x+h).^2 - (x).^2) ./ h;
        case 'central'
            grad = seq .* ((x+h).^2 - (x-h).^2) ./ (2.*h);
        otherwise
            error('Requested difference not available')
    end
end

