function res = SimulateMDScan(ts,feyn,Omega,graph,inhom)

axnum = length(ts);
axlen = zeros(1,axnum);
for i = 1:axnum
    axlen(i) = length(ts{i});
    sh = ones(1,axnum);
    sh(i) = length(ts{i});
    ts{i} = reshape(ts{i},sh);
end

res = zeros(axlen);

if exist('inhom','var')
    sigma = zeros(length(inhom.sigma),axnum);
end

for i = 1:length(feyn)
    feyncurr = feyn{i};
    resn = 1;
    
    for n = 1:length(feyncurr.sign)-1
        ketind = find(matches(graph.Nodes.Name,feyncurr.ket{n+1}));
        braind = find(matches(graph.Nodes.Name,feyncurr.bra{n+1}));

        if exist('inhom','var')
            sigma(:,n) = inhom.sigma .* (inhom.key(:,ketind) - inhom.key(:,braind));
        end
        
        resn = resn * feyncurr.side(n) * 1i .* exp(-1i * Omega(ketind,braind) * ts{n});                          %#ok<FNDSB> 
    end
    
    if exist('inhom','var')
        res = res + resn .* exp(-1/2 * cellprod(sigma' * inhom.rho * sigma));
    else
        res = res + resn;
    end
end

    function ires = cellprod(arr)
        ires = 0;
        for k = 1:axnum
            for l = 1:axnum
                ires = ires + ts{k} .* arr(k,l) .* ts{l};
            end
        end
    end

end