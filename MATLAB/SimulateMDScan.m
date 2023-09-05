function res = SimulateMDScan(ts,feyn,Omega,graph)

axnum = length(ts);
axlen = zeros(1,axnum);
for i = 1:axnum
    axlen(i) = length(ts{i});
    sh = ones(1,axnum);
    sh(i) = length(ts{i});
    ts{i} = reshape(ts{i},sh);
end

res = zeros(axlen);

for i = 1:length(feyn)
    feyncurr = feyn{i};
    resn = 1;
    
    for n = 1:length(feyncurr.sign)-1
        ketind = find(matches(graph.Nodes.Name,feyncurr.ket{n+1}));
        braind = find(matches(graph.Nodes.Name,feyncurr.bra{n+1}));

        resn = resn * feyncurr.side(n) * 1i .* exp(-1i * Omega(ketind,braind) * ts{n});                          %#ok<FNDSB> 
    end
    
    res = res + resn;
end

end