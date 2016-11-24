function varargout = generate_LIF_network_adjacency_matrix(network_params)
addpath ../../../../Research/Sleep_Wake_Matlab_Codes_2013/test_network_data/scratch
addpath ../../../../Research/Sleep_Wake_Matlab_Codes_2013/power-law-estimator/
%%
network_type = network_params{1}; 	% Network type
en = network_params{2}{1};         	% # of excitatory neurons (pos. int.)
in = network_params{2}{2};         	% # of inhibitory neurons (neg. int.)
% ------------------------------------------------------------------------
n = en+(-1)*in;
% ------------------------------------------------------------------------
if strcmp(network_type, 'ER')
    % Edge probability; form edge with probability p
    edge_probability = network_params{3}{1};
    % Consider only the upper triangular section
    adjacency_matrix = triu(rand(n)<edge_probability, 1);
    % Sum with its transpose to guarantee symmetricity
    adjacency_matrix = adjacency_matrix + adjacency_matrix';
elseif strcmp(network_type, 'SF')
    % Power-law exponent alpha
    alpha = network_params{3}{1};
    adjacency_matrix = scale_free_graph_generation(n, alpha, ...
        'zero degree is ok');
else
    display(network_type);
    error('Unexpected inputs!!!');
end
% ------------------------------------------------------------------------
% Average number of neighbors of an excitatory neuron
eneff = sum(sum(adjacency_matrix(:,1:en)))/n;
% Average number of neighbors of an inhibitory neuron
ineff = sum(sum(adjacency_matrix(:,en+1:end)))/n;
% Number of excitatory and inhibiotry neighbors
ei_nbrs = zeros(n,2);
for i = 1:n
    ei_nbrs(i,1) = sum(adjacency_matrix(i, 1:en));
    ei_nbrs(i,2) = sum(adjacency_matrix(i, en+1:n));
end
% ------------------------------------------------------------------------
varargout{1} = adjacency_matrix;
varargout{2} = eneff;
varargout{3} = ineff;
varargout{4} = ei_nbrs;
end