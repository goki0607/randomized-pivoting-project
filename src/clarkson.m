function [x,f,its,flag] = clarkson(Ap, bp, xp, cp)
    global Ag xg bg cg iter max_its;
    Ag = Ap;
    xg = xp;
    bg = bp;
    cg = cp;
    iter = 0;
    max_its = 10*(size(Ap,1)+size(Ap,2));
    S = (1:length(Ap))';
    [x,f,~,flag] = xr(S);
    its = iter;
end

function [x,f,its,flag] = xr(S)
    global Ag bg iter max_its;
    n = length(S);
    d = size(Ag,2);
    Cd = 9*d*d;
    if iter > max_its
        flag = -1;
        return
    elseif n <= Cd
        [x,f,its,flag] = solve(S);
        iter = iter + its;
    else
        Vstar = [];
        while 1
            r = floor(d*sqrt(n));
            R = randsample(setdiff(S,Vstar),r);
            [x,f,its,flag] = xr(union(R,Vstar));
            disp(iter)
            disp(f)
            if flag ~= 0
                return
            end
            V = S(Ag(S,:)*xstar-bg(S,:)<0);
            if length(V) <= 2*sqrt(n)
                Vstar = union(Vstar,V);
            end
            if isempty(V)
                break;
            end
        end
        flag = 0;
    end
end

function [x,f,its,flag] = solve(S)
    global Ag xg bg cg;
    [x,f,its,flag] = evallp(Ag(S,:),bg(S,:),cg,xg,0);
end

function [x,f,its,flag] = evallp(A,b,c,~,option)
% Two phase linear programming solver. This time we wish to solve the
% problem given an A b and c rather than using a mps file. Pption triggers
% what kind of solution strategy we would like to proceed with.
% Inputs:       p -- mps object that describes the problem,
%          option -- which solving approach we would like, between 0 and 3
%                    is the pivoting strategies described in simplex.m, 
%                    4 is random facet and 5 is clarkson's algorithm.
% Outputs:    x -- an (possibly) optimal solution to the given problem,
%             t -- time statistics using timeit,
%             f -- the optimal solution value or -Inf/Inf if unbounded/
%                  infeasible,
%           its -- number of iterations performed by the strategy taken.
%          flag -- exit flag: 2 if unbounded, 1 if infeasible, 0 if 
%                  solution found and -1 if maximum number of iterations.
    p1 = @() phase11(A,b);
    [x,w,flag,i] = p1();
    if flag == 2
        x = (1:d)' + -Inf;
        f = -Inf;
        its = i;
        return
    elseif flag == 1
        its = i;
        x = (1:d)' + Inf;
        f = Inf;
        return
    else
        i1 = i;
    end
    if option >= 0 && option <= 3
        g = @() simplex(w,A,b,x,c,option);
    elseif option == 4
        g = @() simplex_random_facet(w,A,b,x,c);
    elseif option == 5
        error("can't call clarkson within clarkson")
    end
    [x,f,~,its,flag] = g();
    its = its + i1;
end

function [x,w,flag,its] = phase11(A,b)
% Phase 1 LP. Attempts to find a feasible point for the problem defined by
% A (ineqaulity constraints) and b (the right side of the inequality) by
% minimizing the sum of the constraints that are violated. If the problem
% is infeasible the program will have exit flag 1, if the maximum number of
% iterations is reached by the solver the program will have exit flag -1
% and otherwise if a solution is found the program exits with 0. The
% program also returns the working set w to continue onto a possible phase
% 2, the possibly feasible points xnew and the iterations performed by the
% solving strategy. The solving strategy is defined by the option variable
% and described in evallp.m.
    [n,d] = size(A);
    [W0,idx] = licols(A');
    x0 = W0' \ b(idx,:);
    S = A*x0-b;
    V = find(S < 0);
    nv = length(V);
    P = zeros(n,nv);
    P(sub2ind(size(P),V,(1:nv)')) = 1;
    Ap = [A P; zeros(nv,d) eye(nv,nv)];
    bp = [b; zeros(nv,1)];
    cp = sparse([zeros(d,1); ones(nv,1)]);
    [xf,f,e,o] = linprog(cp,-Ap,-bp);
    its = o.iterations;
    if e == 0
        flag = -1;
    elseif f ~=0
       flag = 1;
     else
       flag = 0;
     end
    x = xf(1:d);
    wc = A*x==b;
    [B,~] = licols(A(wc,:)');
    [~,w] = intersect(A,B','rows');
end