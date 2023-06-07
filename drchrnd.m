function r = drchrnd(a,n)
    % take a sample from a dirichlet distribution
    p = length(a);
    r = gamrnd(repmat(reshape(a,1,p),n,1),1,n,p);
    r = r ./ repmat(sum(r,2),1,p);
end