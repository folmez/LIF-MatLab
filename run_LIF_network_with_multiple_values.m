function varargout = run_LIF_network_with_multiple_values(varargin)
%% Run simulation or load from saved file
if strcmp(varargin{1}, 'load from file')
    %% Load from saved file
    load(varargin{2});
else
    %% Network parameters
    type = varargin{1};     % Type of network: 'ER'
    en = varargin{2}(1);    % # of excitatory neurons (positive integer)
    in = varargin{2}(2);    % # of inhibitory neurons (negative integer)
    t_end = varargin{3}(1); % in milliseconds
    dt = varargin{3}(2);    % in milliseconds
    
    network_params_vec = 0.5;   % Other network parameters
    
    input_rate_vec = 4000;  % Poisson input rate
    fS_vec = 0.005;         % Strength of Poisson input

    extI_vec = 0;
    
    tref = 0;               % Refractory period
    nr_sims = 1;
    
    save_workspace = 0;
    plot_stuff = 1;
    % --------------------------------------------------------------------
    i = 4;
    while i<=length(varargin),
        switch varargin{i},
            case 'input_rate_vec'
                input_rate_vec = varargin{i+1};
            case 'input_strength_vec'
                fS_vec = varargin{i+1};
            case 'network_params_vec'
                network_params_vec = varargin{i+1};
            case 'extI_vec'
                extI_vec = varargin{i+1};
            case 'tref'
                tref = varargin{i+1};
            case 'nr_sims'
                nr_sims = varargin{i+1};
            case 'rate'
                if length(varargin{i+1})==1
                    startPR = varargin{i+1}(1);
                    endPR = varargin{i+1}(1);
                    stepPR = 100;
                elseif length(varargin{i+1})==3
                    startPR = varargin{i+1}(1);
                    endPR = varargin{i+1}(2);
                    stepPR = varargin{i+1}(3);
                end
                input_rate_vec = (startPR:stepPR:endPR)';
            case 'save_workspace'
                save_workspace = varargin{i+1};
            case 'plot_stuff'
                plot_stuff = varargin{i+1};
            otherwise,
                display(varargin{i});
                error('Unexpected inputs!!!');
        end
        i = i+2;
    end
    % --------------------------------------------------------------------
    % Organize input vectors
    if length(input_rate_vec)>1
        nr_iters = length(input_rate_vec);
        fS_vec = ones(nr_iters, 1)*fS_vec;
        network_params_vec = ones(nr_iters, 1)*network_params_vec;
        extI_vec = ones(nr_iters, 1)*extI_vec;
        iter_vec = input_rate_vec;
        iter_varname = 'Poisson Rate (#/sec)';
        filename_prefix = 'input_rate_varied';
    elseif length(fS_vec)>1
        nr_iters = length(fS_vec);
        input_rate_vec = ones(nr_iters, 1)*input_rate_vec;
        network_params_vec = ones(nr_iters, 1)*network_params_vec;
        extI_vec = ones(nr_iters, 1)*extI_vec;
        iter_vec = fS_vec;
        iter_varname = 'Poisson strength';
        filename_prefix = 'input_strength_varied';
    elseif length(network_params_vec)>1
        nr_iters = length(network_params_vec);
        fS_vec = ones(nr_iters, 1)*fS_vec;
        input_rate_vec = ones(nr_iters, 1)*input_rate_vec;
        extI_vec = ones(nr_iters, 1)*extI_vec;
        iter_vec = network_params_vec;
        iter_varname = 'Network parameter';
        filename_prefix = 'network_params_varied';
    elseif length(extI_vec)>1
        nr_iters = length(extI_vec);
        fS_vec = ones(nr_iters, 1)*fS_vec;
        input_rate_vec = ones(nr_iters, 1)*input_rate_vec;
        network_params_vec = ones(nr_iters, 1)*network_params_vec;
        iter_vec = extI_vec;
        iter_varname = 'External input';
        filename_prefix = 'external_input_varied';
    else
        warning('Running a single simulation...');
    end
    % --------------------------------------------------------------------
    plot_LIF_network_results = 0;
    display_LIF_network_input_output_summary = 0;
    % --------------------------------------------------------------------
    %% Model
    % Initiate a matrix to record population firing rates:
    % [total, excitatory, inhibitory]
    % Population Firing Rate Matrix
    pfrm = zeros(nr_iters, nr_sims, 3);
    tSIM = tic;
    for i = 1:nr_iters
        for j = 1:nr_sims
            [~, ~, pfrm(i,j,:)] = LIF_network(type, [en in], ...
                network_params_vec(i), [t_end dt], ...
                'external_spike_rate', input_rate_vec(i), ...
                'external_spike_strength', fS_vec(i), ...
                'extI', extI_vec(i), ...
                'tref', tref, ...
                'plot_results', plot_LIF_network_results, ...
                'display_input_output_summary', ...
                display_LIF_network_input_output_summary);
            fprintf('[%3.2fm]: %i/%i completed.\n', ...
                toc(tSIM)/60, (i-1)*nr_sims+j, nr_iters*nr_sims);
        end
    end
    %% Save workspace
    if save_workspace
        filename = ['saved_workspaces/LIF_' filename_prefix ...
            datestr(now,'_mmmdd_HHMMSS') '.mat'];
        save(filename);
    end
end
%%  Plot Results
if plot_stuff
    figure,
    if nr_sims>1
        % Errorbars
        errorbar(iter_vec, mean(pfrm(:,:,2), 2), std(pfrm(:,:,2), 0, 2));
        hold on;
        errorbar(iter_vec, mean(pfrm(:,:,3), 2), std(pfrm(:,:,3), 0, 2));
    else
        % Regular plot
        plot(iter_vec, mean(pfrm(:,:,2), 2));
        hold on;
        plot(iter_vec, mean(pfrm(:,:,3), 2));
    end
    % plot(iter_vec, pop_firing_rate_mat(:,2:3), '.-', 'MarkerSize', 20);
    h_leg = legend('exc.f.r.', 'inh.f.r.', 'Location', 'Best');
    set(h_leg, 'FontSize', 25);
    xlabel(iter_varname);
    ylabel('Population Firing Rate (#/sec)');
end
%% Outputs
varargout{1} = pfrm;
varargout{2} = iter_vec;
varargout{3} = iter_varname;
end