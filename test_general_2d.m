% this script will implement the Lyapunov function learning method on a
% general triangularization of the space, in this code. We consider the 2d
% space
clear all; close all; clc
%% basic configuration
dyn = @(x) [-.9*sin(x(1,:)).*cos(x(2,:))+.2*x(1,:).*x(2,:)+.25*x(2,:).^2;-1.*sin(x(2,:)).*abs(x(1,:)+.2)+.5*x(1,:).*(x(2,:))./(1*cos(x(2,:)-.3*x(1,:)))-.0*x(2,:)];rng(124)

M  = 1.2;  ... (0.4 to 1 choose 1.3, .1 to .4 choose 1.2)

% sample data
num_samples =210;
intv_samp = [-1,1];
intv = [-.4,.4];      ... the region of interest in the stability analysis (along one axis, symmetric)
intv_stb = [-0.1,0.1];      ... the region that is stable
% add some structural nodes to better control the resolution of the grid
resolution = .1;       ... 0.4 to 1 choose .15, .1 to .4 choose .1
[temp_x,temp_y] = meshgrid(intv(1)+resolution:resolution:intv(2)-resolution,intv(1)+resolution:resolution:intv(2)-resolution);
nodes = [temp_x(:),temp_y(:)];
temp = [];
for i = 1:size(nodes,1)
    if nodes(i,1)>=intv_stb(1)-1e-3 && nodes(i,1)<=intv_stb(2)+1e-3 && ...
            nodes(i,2)>=intv_stb(1)-1e-3 && nodes(i,2)<=intv_stb(2)+1e-3
        continue;
    else
        temp(end+1) = i;
    end
end
nodes = nodes(temp,:);

% ------ setup of the gridding ------
num_grid_inner = 5;    ... number of the grid points on the edge of the a-priori basin (0.4 to 1 choose 10, .1 to .4 choose 5)
num_grid_outer = 10;    ... number of the grid points on one edge of the boundary   (0.4 to 1 choose 20, .1 to .4 choose 10)

num_node = 200;     ... number of grid node samples in the region of interest

% grid generation
temp = 0;
while true
    temp_nodes = (intv(2)-intv(1))*rand(1,2)+intv(1);
    if temp_nodes(1)>=intv_stb(1)-1e-2 && temp_nodes(1)<=intv_stb(2)+1e-2 && ...
            temp_nodes(2)>=intv_stb(1)-1e-2 && temp_nodes(2)<=intv_stb(2)+1e-2
        continue;
    else
        nodes(end+1,:) = temp_nodes;
        temp = temp+1;
    end

    if temp>=num_node
        break;
    end
end

idx_int = [1,length(nodes)];

% append the outer boundary to the nodes
temp = linspace(intv(1),intv(2),num_grid_outer)';
temp_nodes = [kron(intv',ones(num_grid_outer,1)),kron(ones(2,1),temp)];
temp_nodes = [temp_nodes;[kron(ones(2,1),temp(2:end-1)),kron(intv',ones(num_grid_outer-2,1))]];
nodes = [nodes;temp_nodes];
% index of the points in the inner boundary
idx_outer = [size(nodes,1)-size(temp_nodes,1)+1,size(nodes,1)];

% append the inner boundary to the nodes
temp = linspace(intv_stb(1),intv_stb(2),num_grid_inner)';
temp_nodes = [kron(intv_stb',ones(num_grid_inner,1)),kron(ones(2,1),temp)];
temp_nodes = [temp_nodes;[kron(ones(2,1),temp(2:end-1)),kron(intv_stb',ones(num_grid_inner-2,1))]];
nodes = [nodes;temp_nodes];
% index of the points in the outer boundary
idx_inner = [size(nodes,1)-size(temp_nodes,1)+1,size(nodes,1)];

% ------ data from the dynamics ------
data = struct;
data.x =  (intv_samp(2)-intv_samp(1)+.1)*(rand(2,num_samples))+intv_samp(1)-.05;
data.dx = dyn(data.x);

figure(1); clf; hold on;
quiver(data.x(1,:),data.x(2,:),data.dx(1,:),data.dx(2,:));
plot(data.x(1,:),data.x(2,:),'o');
axis([intv_samp(1)-0.3, intv_samp(2)+0.3,intv_samp(1)-0.3, intv_samp(2)+0.3])
for i = 1:1:length(data.x)
    [~,y_temp]=ode23tb(@(t,x) dyn(x),[0 5],data.x(:,i));   
    plot(y_temp(:,1),y_temp(:,2));
end

%% ====== step 1: generate the tessellation by triangularization ======
% idxex mapping to the nodes forming the delaunay triangularization
tri_idx = delaunay(nodes);
% remove the triangle in the basin
temp = (idx_inner(1)<=tri_idx).*(tri_idx<=idx_inner(2));
idx = find(sum(temp,2)<3);

tri_idx = tri_idx(idx,:);

% ------ validity of the grid ------
flag = valid_check_2d(tri_idx,nodes,data,M);
if ~flag
    error('The triangulation cannot learn the Lyapunov function!');
end

% visualization of the triangularization
figure(2);clf;hold on
triplot(tri_idx,nodes(:,1),nodes(:,2))

%% ====== step 2: choose the data that is relevant for each sub grids ======
fprintf('select data that are relevant for robust stability constraints...\n')
data_idx = data_refinement(tri_idx,nodes,data,M);

%% ====== step 3: define the optimization problem =======
max_lya = 10;   ... maximal evaluation of the lyapunov function in the outer boundary
min_lya =1;   ... minimal evaluation of the lyapunov function in the inner boundary
mode = 2;     ... lower bound not exact (0.4 to 1 choose 4, .1 to .4 choose 2)
tol = 1e-3;     ... tolerance to enforce strict negative constriant

fprintf('define the Lyapunov function evaluation constraints ...\n')
[cons_eval,g,b,val_vtx] = eval_cons(tri_idx,nodes,idx_inner,idx_outer,idx_int,min_lya,max_lya,tol,mode);
fprintf('define the robust stability constraints ...\n')
[g_aux,t_aux,slack,cons_rb] = rb_stb_relaxed(tri_idx,nodes,data_idx,data,g,M,tol);

fprintf('Finish definition of all the constraints!\n');

obj = sum(sum(slack));

ops =  sdpsettings('solver','mosek','verbose',0);
fprintf('solving the robust optimization problem ...\n')

diagnostics = optimize([cons_eval;cons_rb],obj,ops);
if diagnostics.problem ~= 0 
    error(diagnostics.info);
elseif double(obj) <= 1e-5 - tol*size(tri_idx,1)
    fprintf('Lyapunov function found!\n')
else
    fprintf('Relaxed solution found!\n')
end

% %% plot the slack variables
% figure(3);clf; hold on;
% for i = 1:length(tri_idx)
%     temp = Polyhedron('V',nodes(tri_idx(i,:),:));
%     temp.plot('alpha', 0.4, 'color',[1 0 1]*min(abs(mean(double(slack(i,:))+1e-3))/(1e-3),1));
% end


%% final plot
val_vtx = double(val_vtx);
g = double(g);
b = double(b);
if mode < 3
    % ====== visualization of the learnt function ======
    figure(4);clf; hold on; 
    for i = 1:length(tri_idx)
        temp_vtx = nodes(tri_idx(i,:),:);
        temp_val = double(val_vtx(tri_idx(i,:)));
        temp = Polyhedron('V',[temp_vtx,vec(temp_val)]);
        temp.plot('alpha', 0.4, 'color','b');
    end
    temp_sample = [];
    temp_val_sample = [];
    % ====== sample based validation ======
    for i = 1:length(tri_idx)
        temp_vtx = nodes(tri_idx(i,:),:);
        temp = drchrnd([1,1,1],100);
        temp = temp*temp_vtx;
        temp_sample = [temp_sample;temp];
        temp = double(g(i,:)*dyn(temp'));
        temp_val_sample = [temp_val_sample;temp'];
    end
    figure(5);  clf;hold on; 
    scatter(temp_sample(:,1),temp_sample(:,2),15,temp_val_sample,'filled');
    temp_ind = find(temp_val_sample>=0);
    scatter(temp_sample(temp_ind,1),temp_sample(temp_ind,2),[],zeros(length(temp_ind),3),'filled');
    colorbar;
elseif mode >=3
    figure(4);clf; hold on; 
    % ====== find the largest sublevel set ======
    alpha = min(double(val_vtx(idx_outer(1):idx_outer(2))))-1e-1;
    temp_val_sample = [];
    temp_sample = [];
    temp_poly = {};
    for i = 1:length(tri_idx)
        temp_vtx = nodes(tri_idx(i,:),:);
        vtx = [];
        temp_val = double(val_vtx(tri_idx(i,:)));
        temp_val_vtx = [];
        if min(temp_val)<=alpha-1e-3
            if max(temp_val)>alpha
                [temp_min_val,temp_idx] = sort(temp_val);
                if temp_min_val(2)<=alpha
                    temp_idx = temp_idx(1:2);
                    temp_min_val = temp_min_val(1:2);
                else
                    temp_idx = temp_idx(1);
                    temp_min_val = temp_min_val(1);
                end
                for j = 1:3
                    if temp_val(j) <=alpha
                        vtx = [vtx;temp_vtx(j,:)];
                        temp_val_vtx = [temp_val_vtx;temp_val(j)];
                    else
                        for k = 1:length(temp_idx)
                            vtx = [vtx;(alpha-temp_min_val(k))/(temp_val(j)-temp_min_val(k))*(temp_vtx(j,:)-temp_vtx(temp_idx(k),:))+temp_vtx(temp_idx(k),:)];
                            temp_val_vtx = [temp_val_vtx; alpha];
                        end
                    end
                end
            else
                vtx = temp_vtx;
                temp_val_vtx = temp_val;
            end

            temp = drchrnd(ones(1,size(vtx,1)),1500);
            temp = temp*vtx;
            temp_sample = [temp_sample;temp];
            temp = double(g(i,:)*dyn(temp'));
            temp_val_sample = [temp_val_sample;temp'];
            temp_poly{end+1} = Polyhedron('V',[vtx]);
            temp = Polyhedron('V',[vtx,vec(temp_val_vtx)]);
            temp.plot('alpha', 0.4, 'color','b');
        else
            continue;
        end
    end
    figure(5);clf; hold on; 
    scatter(temp_sample(:,1),temp_sample(:,2),15,temp_val_sample,'filled');
    temp_ind = find(temp_val_sample>=0);
    scatter(temp_sample(temp_ind,1),temp_sample(temp_ind,2),[],zeros(length(temp_ind),3),'filled');
    colorbar;
    for i = 1:length(temp_poly)
        temp_poly{i}.plot('alpha', 0.1, 'color','b');
    end
end


