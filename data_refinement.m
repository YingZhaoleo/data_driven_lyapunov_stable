function data_idx = data_refinement(idx,nodes,data,M)
% for each simplex, only keep the data that is relevant to its analysis of
% stability
for nsplx = 1:size(idx,1)
    v = nodes(idx(nsplx,:),:); ... vertices of this simplex
    temp_ind = [];
    for ndata = 1:size(data.x,2)
            if (norm(data.dx(:,ndata),2)-M*norm(data.x(:,ndata)-v(1,:)',2))>=0 || ...
               (norm(data.dx(:,ndata),2)-M*norm(data.x(:,ndata)-v(2,:)',2))>=0 || ...
               (norm(data.dx(:,ndata),2)-M*norm(data.x(:,ndata)-v(3,:)',2))>=0
                temp_ind(end+1) = ndata;
            end
    end
    data_idx{nsplx} = temp_ind;
        
end

end