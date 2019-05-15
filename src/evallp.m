function [x,t,f,its,flag] = evallp(p,option)
% Two phase linear programming solver. P is a mps object which describes
% the problem we wish to solve while option triggers what kind of solution
% strategy we would like to proceed with.
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
    %p = mpsread(str);
    f = p.f;
    Aineq = -p.Aineq;
    bineq = -p.bineq;
    Aeq = [-p.Aeq; p.Aeq];
    beq = [-p.beq; p.beq];
    d = length(f);
    idx1 = find(p.lb > -Inf);
    idx2 = find(p.ub < Inf);
    I = eye(d,d);
    Alub = [I(idx1,:); -I(idx2,:)];
    blub = [p.lb(idx1,:); -p.ub(idx2,:)];
    A = [Aineq; Aeq; Alub];
    b = [bineq; beq; blub];
    c = f;
    p1 = @() phase11(A,b);
    tic
    [x,w,flag,i] = p1();
    disp('phase 1 finished')
    t = toc;
    if flag == 2
        x = (1:d)' + -Inf;
        f = -Inf;
        %t = timeit(p1);
        its = i;
        return
    elseif flag == 1
        x = (1:d)' + Inf;
        f = Inf;
        its = i;
        %t = timeit(p1);
        return
    else
        %t1 = t;
        i1 = i;
    end
    if option >= 0 && option <= 3
        g = @() simplex(w,A,b,x,c,option);
    elseif option == 4
        g = @() simplex_random_facet(w,A,b,x,c);
    elseif option == 5
        g = @() clarkson(A,b,x,c);
    end
    if option == 5
        tic
        [x,f,its,flag] = g();
        t1 = toc;
    else
        tic
        [x,f,~,its,flag] = g();
        t1 = toc;
    end
    disp('phase 2 finished')
    %t2 = timeit(g);
    t = t1;
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
    %F = find(S >= 0);
    V = find(S < 0);
    %nf = length(F);
    nv = length(V);
    P = zeros(n,nv);
    P(sub2ind(size(P),V,(1:nv)')) = 1;
    Ap = [A P; zeros(nv,d) eye(nv,nv)];
    %Ap2 = [A(F,:) zeros(nf,nv); A(V,:) eye(nv,nv); zeros(nv,d) eye(nv,nv)];
    bp = [b; zeros(nv,1)];
    %bp2 = [b(F,:); b(V,:); zeros(nv,1)];
    cp = sparse([zeros(d,1); ones(nv,1)]);
    %rp = A(V,:)*x0-b(V,:);
    %xp = [x0; -rp];
    %wp = [idx; V];
    [xf,f,e,o] = linprog(cp,-Ap,-bp);
    its = o.iterations;
    if e == 0
        flag = -1;
    elseif f ~=0
       flag = 1;
     else
       flag = 0;
     end
    %{
    if option >= 0 && option <= 3
        [xf,f,w,its,flag] = simplex(wp,Ap,bp,xp,cp,option);
    elseif option == 4
        [xf,f,w,its,flag] = simplex_random_facet(wp,Ap,bp,xp,cp);
    elseif option == 5
        [xf,f,its,flag] = clarkson(Ap,bp,xp,cp);
        w = find(A*xf(1:d)==b);
    end
    %}
    x = xf(1:d);
    %w = w(1:d);
    %[~,idx2] = licols(A(w,:)');
    wc = A*x==b;
    [B,~] = licols(A(wc,:)');
    [~,w] = intersect(A,B','rows');
end