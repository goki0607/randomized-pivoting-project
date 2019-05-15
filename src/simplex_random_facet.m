function [x,f,w,its,flag] = simplex_random_facet(wp, Ap, bp, xp, cp)
% Simplex random facet function that performs the simplex method using 
% the random facet pivoting rule. The maximum number of iterations is
% defined to be the maximum number of pivot moves made after the first
% recursive call and is given by 10*(#ineq+#vars).
% working set is increased.
% Inputs:     wp -- working set at x, must be vector of indices,
%             Ap -- set of constraints, must be a matrix of constraints,
%             bp -- set of inequalities corresponding to set of
%                   constraints, assumed to be in all-inequality form
%                   (... >= b_i),
%             xp -- the initial starting point, must be a vector and a
%                   vertex of the problem,
%             cp -- coefficients on the objective function, must be a
%                   vector, i.e. min c^Tx,
% Outputs:    x -- an optimal solution to the given problem,
%             f -- the optimal solution value or -Inf if unbounded,
%             w -- the working set upon exit of the program,
%           its -- number of iterations performed,
%          flag -- exit flag: 1 if unbounded, 0 if solution found and -1 if
%                  maximum number of iterations reached and solution is not
%                  optimal.
    global wg Ag bg cg iters max_its;
    wg = wp;
    Ag = Ap;
    bg = bp;
    cg = cp;
    iters = 0;
    max_its = 10*(size(Ap,1)+size(Ap,2));
    nums = (1:length(Ap))';
    [w,x,flag] = random_facet(nums,wp,xp);
    its = iters;
    if flag == 1
        f = Inf;
    else
        f = cp'*x;
    end
    w = sort(w);
end

function [B,x,flag] = random_facet(F,B0,x0)
% The random facet algorithm.
    global iters max_its ;
    if iters > max_its
        B = B0;
        x = x0;
        flag = -1;
    elseif isempty(intersect(F,B0))
        B = B0;
        x = x0;
        flag = 0;
    else
        f = random_f(intersect(F,B0));
        [B1,x1,flag1] = random_facet(setdiff(F,f),B0,x0);
        [B2,x2,flag2,optimal] = pivot(F,B1,f,x1);
        if flag1 == 1 || flag2 == 1
            flag = 1;
            B = B0;
            x = x0;
        elseif flag1 == -1 || flag2 == -1
            flag = -1;
            B = B0;
            x = x0;
        elseif optimal == 1
            B = B1;
            flag = 0;
            x = x1;
        elseif ~isequal(B1,B2)
            [B,x,flag] = random_facet(setdiff(F,f),B2,x2);
        else
            B = B1;
            x = x1;
            flag = 0;
        end
    end
end

function f = random_f(F)
% Auxilarry function to obtain a random facet from a set of facets F.
    idx = randi([1 length(F)], 1);
    f = F(idx);
end

function [B,x,flag,optimal] = pivot(~,B1,f,x0)
% Function to pivot away from constraint f. We randomly pivot to another
% vertex.
    global Ag bg cg iters;
    optimal = 0;
    iters = iters + 1;
    W = Ag(B1,:);
    lambda = W' \ cg;
    s = find(B1==f);
    if (lambda(s) >= -1e-7)
        B = B1;
        x = x0;
        flag = 0;
        optimal = 1;
        return
    end
    es = zeros(length(lambda),1);
    es(s) = 1;
    p = W \ es;
    [Ap,idx] = setdiff(Ag,W,'rows');
    D = idx(Ap*p<-1e-7);
    if isempty(D)
        B = B1;
        x = x0;
        flag = 1;
        return
    end
    gamma = (Ag(D,:)*x0-bg(D)) ./ (-Ag(D,:)*p);
    alpha_f = min(gamma);
    x = x0 + alpha_f*p;
    S = D(gamma==alpha_f,:);
    if isempty(S)
        B = B1;
        x = x0;
        flag = 1;
        return
    end
    % t = min(S);
    t = random_in(S);
    if t == -1
        B = B1;
        x = x0;
        flag = 1;
        return
    end
    B = B1;
    B(s) = t;
    flag = 0;
end

function t = random_in(S)
% Random edge rule for incoming constraint.
    if isempty(S)
        t = -1;
        return
    end
    idx = randi([1 size(S,1)], 1);
    t = S(idx);
end