function [cons_eval,g,b,val_vtx] = eval_cons(idx,nodes,idx_inner,idx_outer,idx_int, min_lya,max_lya,tol,mode)
% define the evaluation constraint on the Lyapunov learning problem,
% including the continuity and the sub level set, this code works for 2d
% Arguments:
%     idx_inner/outer/int: the starting/tail index of verteces indeces that lie on the boundary of the 
%             prior forward stable set Xs, the boudnary of set to be analyzed X, and the interior of X\X_s
%             Meanwhile, idx is n x 3 matrix, each row is the index for one specific simplex
%     mode: 1: lower bound of the lypapunov function is exact,upper bound of the lyapunov function is exact
%         2: lower bound of the lypapunov function less than min_lya,upper bound of the lyapunov function is exact
%         3: lower bound of the lypapunov function less than min_lya,upper bound of the lyapunov function larger than max_lya
% Return:
%     g,b: the parameters defines each affine pieces
%     val_vtx: the evaluation on each vertex
%     cons_eval: the yalmip constraint object including all these constraints

g = sdpvar(size(idx,1),2,'full');
b = sdpvar(size(idx,1),1,'full');
val_vtx = sdpvar(size(nodes,1),1,'full');

% some constraint on the inner boundary to improve the regularity
if mode ==1
    % forced to be constant on the inner boundary, used in the sequential
    % algorithm
    cons_eval = [val_vtx(idx_inner(1):idx_inner(2))== min_lya]; 
    cons_eval = [cons_eval;val_vtx(idx_outer(1):idx_outer(2))==max_lya];
elseif mode ==2
    % no need to be constant on the inner boundary
    cons_eval = [val_vtx(idx_inner(1):idx_inner(2))<= min_lya];
    cons_eval = [cons_eval;val_vtx(idx_outer(1):idx_outer(2))==max_lya];
elseif mode ==3
    % no need to be constant on the inner boundary
    cons_eval = [val_vtx(idx_inner(1):idx_inner(2))<= min_lya];
    cons_eval = [cons_eval;val_vtx(idx_outer(1):idx_outer(2))>=max_lya];
elseif mode ==4
    % no need to be constant on the inner boundary
    cons_eval = [val_vtx(idx_inner(1):idx_inner(2))<= min_lya];
    cons_eval = [cons_eval;val_vtx(idx_outer(1):idx_outer(2))>=min_lya+tol];
else
    error('Unknonw boundary condition for function evaluation.')
end

% sub-level set constraint on the outer boundary
cons_eval = [cons_eval;val_vtx(idx_int(1):idx_int(2))<=max_lya-tol];

for nsplx = 1:size(idx,1)
    v = nodes(idx(nsplx,:),:); ... vertices of this simplex

    % unroll the evaluation cons on all the vertices
    cons_eval = [cons_eval;val_vtx(idx(nsplx,1))==g(nsplx,:)*v(1,:)'+b(nsplx)];
    cons_eval = [cons_eval;val_vtx(idx(nsplx,2))==g(nsplx,:)*v(2,:)'+b(nsplx)];
    cons_eval = [cons_eval;val_vtx(idx(nsplx,3))==g(nsplx,:)*v(3,:)'+b(nsplx)];
end

end