%% Dec 13, 2023 - Kelsey Bates
% Find spectra corresponding to a set of Feynman diagrams
%
% Input:
%   ts:    A cell array of the time delays between each pulse, in order.
%          Each interaction should be an element in the array, so 
%          length(ts) = length(feyn{1}.sign), for instance. Time delays
%          which are varied can be an array of times.
%   feyn:  Cell array of Feynman diagrams, output of FeynmanFinderFunc.
%   Omega: A term corresponding to the evolution of the density matrix. Can
%          include both the oscillation frequency omega and the dephasing
%          parameter i * gamma.
%   graph: Directed graph of energy levels. Should be the same as the one
%          used to generate feyn. Only used to get indices for the bras and
%          kets, could remove in the future.
%   inhom: Optional argument for inhomogeneous broadening. A struct with
%          three fields. These fields will refer to the base energies,
%          which are the energy differences which make up the system, with
%          some examples given below.
%            rho:    Correlation matrix. Each row and column corresponds to
%                    one of the base energies. Matrix should have ones on
%                    the diagonal. Off-diagonal elements represent the
%                    degree to which those energies are correlated.
%            sigma:  An array of the inhomogeneous linewidths of the base
%                    energies. Note that if we turn sigma into a diagonal
%                    matrix, then sigma * rho * sigma should be positive
%                    definite.
%            key:    A matrix relating the base energies with the states.
%          Example 1:  Consider a three-level V system with one ground
%                      state g and two excited states e1 and e2. Here there
%                      are two base energies, E_e1-E_g and E_e2-E_g. For
%                      this system our key will be [0 1 0; 0 0 1]:
%                                   g e1 e2
%                        E_e1-E_g  [0  1  0;
%                        E_e2-E_g   0  0  1]
%          Example 2:  Consider a diamond system with one ground state gg,
%                      two singly excited states ge and eg, and one doubly
%                      excited state ee. Here, there are two base energies,
%                      E_ge-E_gg~=E_ee-E_eg and E_eg-E_gg~=E_ee_E_ge. For
%                      this system our key will be [0 1 0 1; 0 0 1 1]:
%                                    gg ge eg ee
%                        E_ge-E_gg  [ 0  1  0  1;
%                        E_eg-E_gg    0  0  1  1]
%
% Output:
%   res:  A complex N-dimensional array, where N is the length of ts. Each
%         point corresponds to the value at a different collection of time
%         delays. The function squeeze(res) can be used to remove
%         dimensions of length 1.

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
    resn = ones(axlen);
    
    for n = 1:length(feyncurr.sign)-1
        ketind = find(matches(graph.Nodes.Name,feyncurr.ket{n+1}));
        braind = find(matches(graph.Nodes.Name,feyncurr.bra{n+1}));

        if exist('inhom','var')
            sigma(:,n) = inhom.sigma .* (inhom.key(:,ketind) - inhom.key(:,braind));
        end
        
        resn = resn * feyncurr.side(n) * 1i .* exp(-1i * Omega(ketind,braind) * ts{n});
    end
    
    if exist('inhom','var')
        res = res + 1i * feyncurr.side(end) * resn .* exp(-1/2 * cellprod(sigma' * inhom.rho * sigma));
    else
        res = res + 1i * feyncurr.side(end) * resn;
    end
end

    function ires = cellprod(Sigma)
        ires = zeros(axlen);
        for k = 1:axnum
            for l = 1:axnum
                ires = ires + ts{k} .* Sigma(k,l) .* ts{l};
            end
        end
    end

end
