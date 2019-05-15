function [x,f,w,its,flag] = simplex(w, A, b, x, c, pivot)
% Simplex method function that performs the simplex method using a pivoting
% code which is an integer between 0 and 3. 0 is the textbook rule, 1 is
% Bland's rule, 2 is steepest edge and 3 is random edge. The maximum number
% of iteratons is determined by the Matlab linprog rule 10*(# of equality
% constraints + # of inequality constraints + # of variables) which in our
% case becomes 10*(#ineq+#vars).
% working set is increased.
% Inputs:      w -- working set at x, must be vector of indices,
%              A -- set of constraints, must be a matrix of constraints,
%              b -- set of inequalities corresponding to set of
%                   constraints, assumed to be in all-inequality form
%                   (... >= b_i),
%              x -- the initial starting point, must be a vector and a
%                   vertex of the problem,
%              c -- coefficients on the objective function, must be a
%                   vector, i.e. min c^Tx,
%          pivot -- must be in the range [0,3] as explained.
% Outputs:    x -- an optimal solution to the given problem,
%             f -- the optimal solution value or -Inf if unbounded,
%             w -- the working set upon exit of the program,
%           its -- number of iterations performed,
%          flag -- exit flag: 1 if unbounded, 0 if solution found and -1 if
%                  maximum number of iterations reached and solution is not
%                  optimal.
	format compact
	format long e
	W = A(w,:);
	lambda = W' \ c;
	k = 1;
    max_its = 10*(size(A,1)+size(A,2));
	while any(lambda<-1e-7) && k <= max_its  %#ok<ALIGN>
        if pivot==0
            s = bland_out(w, lambda);
        elseif pivot==1
            s = dantzig_out(lambda);
        elseif pivot==2
            s = steepest_out(lambda, W);
        elseif pivot==3
            s = random_out(lambda);
        else
            error('please give a correct pivot code')
        end
        es = zeros(length(lambda),1);
        es(s) = 1;
        p = W \ es;
        [Ap,idx] = setdiff(A,W,'rows');
        D = idx(Ap*p<-1e-7);
        if isempty(D)
            flag = 1;
            f = -Inf;
            its = k;
            return
        end
        gamma = (A(D,:)*x-b(D)) ./ (-A(D,:)*p);
        alpha_f = min(gamma);
        alpha = alpha_f;
        x = x + alpha*p;
        S = D(gamma==alpha_f,:);
        if pivot==0
            t = bland_in(S);
            if t == -1
                flag = 1;
                f = -Inf;
                its = k;
                return
            end
        elseif pivot==1
            t = dantzig_in(A, S, p);
        elseif pivot==2
            t = steepest_in(A, S, p);
        elseif pivot==3
            t = random_in(S);
            if t == -1
                flag = 1;
                f = -Inf;
                its = k;
                return
            end
        else
            error('please give a correct pivot code')
        end
        w(s) = t;
        W = A(w,:);
        lambda = W' \ c;
        k = k + 1;
    end
    if k > max_its
        if any(lambda<-1e-7)
            flag = -1;
        else
            flag = 0;
        end
    else
        flag = 0;
    end
    f = c' * x;
    its = k;
    w = sort(w);
end

function s = bland_out(w, lambda)
% Bland's rule for outgoing constraint.
    choices = lambda<0;
    elem = min(w(choices));
    s = find(w==elem);
end

function t = bland_in(S)
% Bland's rule for incoming constraint.
    if isempty(S) 
        t = -1;
        return
    end
    t = min(S);
end
    
function s = dantzig_out(lambda)
% Danzig's rule for outgoing constraint.
    [~, s] = min(lambda);
end

function t = dantzig_in(A, S, p)
% Danzig's rule for incoming constraint.
    [~, idx2] = min(A(S,:)*p);
    t = S(idx2);
end

function s = steepest_out(lambda, W)
% Steepest edge rule for outgoing constraint.
    idx = find(lambda<0);
    lambdas = zeros(length(idx),1);
    for i=1:length(idx)
        es = zeros(length(lambda),1);
        es(idx(i)) = 1;
        p = W \ es;
        lambdas(i) = lambda(i)/norm(p);
    end
    [~,idx2] = min(lambdas);
    s = idx(idx2);
end

function t = steepest_in(A, S, p)
% Steepest edge rule for incoming constraint (nothing fancy).
    [~, idx2] = min(A(S,:)*p);
    t = S(idx2);
end

function s = random_out(lambda)
% Random edge rule for outgoing constraint.
    choices = find(lambda<0);
    idx = randi([1 size(choices,1)], 1);
    s = choices(idx);
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