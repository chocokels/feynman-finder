%% Sept 12, 2019 - Kelsey Bates
% Recursive function to find Feynman diagrams of a system.
%
% Input:
%   ket:        Name of starting ket
%   bra:        Name of starting bra
%   graph:      Directed graph of energy levels. Directions should point
%               from ground to excited states. The names of the nodes
%               (graph.Nodes.Name) should include the names of the starting
%               bra and ket.
%   eta:        Vector array of signs of signal pathways. E.g. [-1 1 1 -1]
%               for S1.
%   heterodyne: 0 for photoluminescence, 1 for heterodyne. Note that
%               heterodyne requires last eta to be -1 to work properly.
%
% Output:
%   feyn: Cell array of structs. Each index in the cell array corresponds
%         to a single Feynman diagram. The struct has four fields.
%           ket:  Cell array of the kets in the Feynman diagram (left side
%                 states).
%           bra:  Cell array of the bras in the Feynman diagram (right side
%                 states).
%           sign: A vector array of the signs of the signal pathway
%                 (essentially just a copy of eta).
%           side: The side of the interaction on the Feynman diagram. 1 is
%                 left and -1 is right.
%         As an example (if you are unfamiliar with cell arrays and
%         structs), to get the second Feynman diagram, use
%           feyn{2}
%         To get the third ket of the second Feynman diagram, use
%           feyn{2}.ket{3}
%         To get the first sign of the second Feynman diagram, use
%           feyn{2}.sign(1)

function feyn = FeynmanFinderFunc(ket,bra,graph,eta,heterodyne)

curr.ket = ket;
curr.bra = bra;
feyn = {};

if isempty(eta)
    if strcmp(ket,bra) && (heterodyne || ~isempty(predecessors(graph,ket)))
        curr.ket = {ket};
        curr.bra = {bra};
        curr.sign = [];
        curr.side = [];
        feyn = {curr};
    end
else
    sign = eta(1);
    side = 1;
    if eta(1) == 1
        states = successors(graph,ket);
    else
        states = predecessors(graph,ket);
    end
    if ~isempty(states)
        for state = 1:length(states)
            next = FeynmanFinderFunc(states{state},bra,graph,eta(2:end),heterodyne);
            if ~isempty(next)
                for i = 1:length(next)
                    next{i}.ket = [ket, next{i}.ket];
                    next{i}.bra = [bra, next{i}.bra];
                    next{i}.sign = [sign, next{i}.sign];
                    next{i}.side = [side, next{i}.side];
                end
                feyn = [feyn;next];
            end
        end
    end
    side = -1;
    if eta(1) == -1
        states = successors(graph,bra);
    else
        states = predecessors(graph,bra);
    end
    if ~isempty(states) && (~heterodyne || (size(eta,2) > 1))
        for state = 1:length(states)
            next = FeynmanFinderFunc(ket,states{state},graph,eta(2:end),heterodyne);
            if ~isempty(next)
                for i = 1:length(next)
                    next{i}.ket = [ket, next{i}.ket];
                    next{i}.bra = [bra, next{i}.bra];
                    next{i}.sign = [sign, next{i}.sign];
                    next{i}.side = [side, next{i}.side];
                end
                feyn = [feyn;next];
            end
        end
    end
end

end