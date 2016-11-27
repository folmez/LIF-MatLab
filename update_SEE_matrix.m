function varargout = update_SEE_matrix(varargin)
%% Input arguments
spike_timing_dependent_plasticity = varargin{1};    % 'on' or 'off'
SEE_matrix = varargin{2};
idx = varargin{3};      % Index set of neurons which fires at the current 
                        % time timestep [t, t+dt]
en = varargin{4};       % # of excitatory neurons
A = varargin{5};        % Adjacency matrix
lft_prev = varargin{6};  % Last firing times of each neuron before t
lft_next = varargin{7};  % Last firing times of each neuron before t+dt
track_synaptic_strength = varargin{8};
synaptic_strength_lb = varargin{9}(1);
synaptic_strength_ub = varargin{9}(2);

display_progress_for_two = track_synaptic_strength;
SE1E2_track = [];

%% Model
if strcmp(spike_timing_dependent_plasticity, 'on')
    % We only consider E-E synapses, ignore inhibitory neuron firings
    e_idx = idx(idx<en+1);
    
    % Order e_idx so that the excitatory neuron that fires first is 
    % considered first in the update loop
    [~, I] = sort(lft_next(e_idx));
    e_idx = e_idx(I);
    
    % Find # of excitatory neuron firings
    nr_exc_firings = length(e_idx);
    
    % Track the syanpses E1->E2 and E2->E1
    track_count = 0;
    SE1E2_track = zeros(3,4); % time-id-e1e2-e2e1
    
    % Update loop
    for i = 1:nr_exc_firings
        % Find the current neuron id
        id = e_idx(i);
        
        % Find its exctitatory neighbors
        id_nbrs = find(A(id, 1:en));
        
        % If neighbors did not fire yet, they must be excluded
        id_nbrs(lft_prev(id_nbrs)==0) = [];
      
        if ~isempty(id_nbrs)
            % Neuron id just fired:
            % (i)  strengths of synapses from its neighbors must be 
            %      potentiated
            SEE_matrix(id_nbrs, id) = SEE_matrix(id_nbrs, id) + ...
                delta_synaptic_strength(SEE_matrix(id_nbrs, id), ...
                lft_next(id)-lft_prev(id_nbrs));
            % (ii) strengths of synapses to its neighbors must be 
            %      depressed
            SEE_matrix(id, id_nbrs) = SEE_matrix(id, id_nbrs) + ...
                delta_synaptic_strength(SEE_matrix(id, id_nbrs), ...
                (lft_prev(id_nbrs)-lft_next(id))');
        end
       
        % Correct strength values that are out of bounds
        SEE_matrix(SEE_matrix<synaptic_strength_lb) = synaptic_strength_lb;
        SEE_matrix(SEE_matrix>synaptic_strength_ub) = synaptic_strength_ub;
                
        % Track and display first two neurons
        if (id==1 && ismember(2, id_nbrs)) || ...
                (id==2 && ismember(1, id_nbrs))
            % Display progress
            if display_progress_for_two
                
                fprintf('\t%1.1f', lft_next(id));
                fprintf('\t%i', id);
                fprintf('\t%1.2e', SEE_matrix(1, 2));
                fprintf('\t%1.2e\n', SEE_matrix(2, 1));
            end
            % Update tracker
            SE1E2_track(track_count+1, :) = ...
                [lft_next(id) id SEE_matrix(1, 2) SEE_matrix(2, 1)];
            track_count = track_count+1;
        elseif ismember(id, [1 2]) && lft_prev(id) == 0
            % Display header
            if display_progress_for_two
                fprintf('  Firing Time');
                fprintf('\tID');
                fprintf('\tE1->E2');
                fprintf('\t\tE2->E1\n');
                fprintf('START:\t%1.1f', lft_next(id));
                fprintf('\t%i', id);
                fprintf('\t%1.2e', SEE_matrix(1, 2));
                fprintf('\t%1.2e\n', SEE_matrix(2, 1));
            end
            % Update tracker
            SE1E2_track(track_count+1, :) = ...
                [lft_next(id) id SEE_matrix(1, 2) SEE_matrix(2, 1)];
            track_count = track_count+1;
        end
        
        % Update last firing time of neuron id
        lft_prev(id) = lft_next(id);
    end
    
    % Remove unused rows from tracker
    SE1E2_track(track_count+1:end,:) = [];
end

%% Output arguments
varargout{1} = SEE_matrix;
varargout{2} = SE1E2_track;
end