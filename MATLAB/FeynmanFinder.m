%% Example energy levels

graph = digraph([1,1],[2,3]);       % Define graph
graph.Nodes.Name = {'0' 'i' 'j'}';  % Name nodes

% graph = digraph([1,2],[2,3]);
% graph.Nodes.Name = {'0' '1' '2'}';

% graph = digraph([1,1,2,3],[2,3,4,4]);
% graph.Nodes.Name = {'0' 'i' 'j' '2'}';

% load('SiVGraph.mat');
% graph = SiVGraph;

%% Plot energy levels
figure(1);clf;
plot(graph);
set(gca,'YDir','reverse');

%% Set eta 
eta = [-1 1 1 -1]; % S1
% eta = [1 -1 1 -1]; % S2
% eta = [1 1 -1 -1]; % S3

%% Find diagrams
feyn = FeynmanFinderFunc('0','0',graph,eta,0);
% feyn = FeynmanFinderFunc('gl','gl',graph,eta,0);  % Use for SiV

%% Print diagrams
PrintFeynman(feyn); 