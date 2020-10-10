%% Sept 12, 2019 - Kelsey Bates
% Prints Feynman diagrams in a human-readable form. (Not identical to a
% hand-drawn diagram, as vectors pointing away from the diagram will be a
% row off. I optimized simplicity over perfection.)
%
% Input:
%   feyn: Output of FeynmanFinderFunc. Note that it requires a cell array
%         input, so use smooth parentheses to index if needed. E.g.:
%           PrintFeynman(feyn);
%           PrintFeynman(feyn(2));
%           PrintFeynman(feyn([2 4]));

function num = PrintFeynman(feyn)

num = length(feyn);
len = length(feyn{1}.sign);
for i = 1:num
    curr = feyn{i};
    disp([' |' curr.ket{end} '><' curr.bra{end} '| ']);
    for j = len:-1:1
        if curr.sign(j) == 1
            symb = '/';
        elseif curr.sign(j) == -1
            symb = '\';
        end
        left = ' ';
        right = ' ';
        if curr.side(j) == 1
            left = symb;
        elseif curr.side(j) == -1
            right = symb;
        end
        disp([left '|' curr.ket{j} '><' curr.bra{j} '|' right]);
    end
    disp(' ');
end

end