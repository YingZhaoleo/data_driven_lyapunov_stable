function [g_aux,t_aux,slack,cons_rb] = rb_stb_relaxed(idx,nodes,data_idx,data,g,M,tol)
% impose the strict Lyapunov decrease constraint
% Return:
%   aux: aux is the auxilary variable to define the SOC constraint
%   slack: the slack variables to define the objective

slack = sdpvar(size(idx,1),3,'full');
cons_rb = [slack(:)>=-tol];

for nsplx = 1:size(idx,1)
    temp_size = length(data_idx{nsplx});
    g_aux{nsplx} = sdpvar(temp_size,2,'full');
    t_aux{nsplx} = sdpvar(temp_size,3,'full');
    cons_rb = [cons_rb;sum(g_aux{nsplx},1)==g(nsplx,:)];

    v = nodes(idx(nsplx,:),:); ... vertices of this simplex
    val = sum(sum(g_aux{nsplx}.*data.dx(:,data_idx{nsplx})')); ... linear part in the robust constraint
    for i = 1:3
        temp_dist = sum((repmat(v(i,:)',1,temp_size)-data.x(:,data_idx{nsplx})).^2,1);
        temp_dist = sqrt(temp_dist);
        for ndata = 1:temp_size
            cons_rb = [cons_rb;cone(M*temp_dist(ndata)*g_aux{nsplx}(ndata,:),t_aux{nsplx}(ndata,i))];
        end
    end
    val = val+sum(t_aux{nsplx},1);
    cons_rb = [cons_rb;slack(nsplx,:)>=val];

end

end