function [Ap,idx] = licols(A,tol)
% Function to retrieve the linearly independent columns of a given matrix A
% along with the set of indices the selected columns with respect to the
% original matrix A. The function uses the QR decomposition of A to find
% all such columns.
% Inputs:    A -- the matrix which we wish to find linearly independent
%                columns for,
%          tol -- the tolerance level for the selections, default is 1e-10.
% Outputs:  Ap -- the matrix which only contains linearly independent
%                columns where the columns of Ap are a subset of the 
%                columns of A,
%          idx -- list of indices of A which are linearly independent.
     if ~nnz(A) %% A is the zero matrix, cannot have independent columns
         Ap = [];
         idx = [];
     else 
        if nargin<2, tol=1e-10; end
        [~,R,E] = qr(A,0); 
        if ~isvector(R)
            dr = abs(diag(R));
        else
            dr = R(1);   
        end
        r = find(dr >= tol*dr(1),1,'last'); %% estimate the rank
        %idx = sort(E(1:r));
        idx = E(1:r);
        Ap = A(:,idx);
        idx = idx';
     end
end