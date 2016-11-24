function test_LIF_network(varargin)
test_type_num = varargin{1};
if length(varargin)>1
    save_workspace = varargin{2};
end

if test_type_num == 4
    %% HALF WITH EXTERNAL INPUTS, OTHER HALF WITHOUT
    %  but with time-dependent external input
    % David's message:
    % Similar to the third case but with a time-dependent external input
    % which is like this: gradual increase from Imin (0) to Imax (100) and
    % slowly decreasing back
    
    if strcmp(varargin{2}, 'load from file')
        load(varargin{3});
        save_workspace = 0;
    else
        % Input parameters
        t_unit = 1000;
        en = varargin{3}(1);
        in = varargin{3}(2);
        edge_prob = 1;
        t_end = varargin{4}(1);
        dt = varargin{4}(2);
        poisson_rate = varargin{5}(1);
        poisson_strength = varargin{5}(2);
        SE = varargin{5}(3);
        SI = varargin{5}(4);
        
        extI_type = 'time-dependent #1';
        Imin = varargin{6}(1);
        Imax = varargin{6}(2);
        nr_extI_steps = varargin{6}(3); % must be even
        
        plot_LIF_results = 1;

        [extI_matrix, extI_step_ends, extI_step_values] = ...
            get_external_input_current_matrix(...
            extI_type, [en in], t_end, dt, [Imin Imax nr_extI_steps]);
        
        % Model
        [~, ~, ~, rast] = LIF_network(...
            'ER', [en in], edge_prob, [t_end dt], 'SE', SE, ...
            'plot_results', plot_LIF_results, ...
            'extI_matrix', extI_matrix, ...
            'rate_strengths', [poisson_rate poisson_strength SE SI]);
        
        % Crunch
        nr_steps = length(extI_step_ends)-1;
        exc_fir_rate_first_half_vec = zeros(nr_steps, 1);
        first_half_idx = rast(:,1)<en/2+1;
        exc_fir_rate_second_half_vec = zeros(nr_steps, 1);
        second_half_idx = rast(:,1)>en/2 & rast(:,1)<=en;
        if in ~= 0
            inh_fir_rate_vec = zeros(nr_steps, 1);
            inh_idx = rast(:,1)>en;
        end
        for i=1:nr_steps
            local_time_window_length = extI_step_ends(i+1) - ...
                extI_step_ends(i);
            local_fir_time_idx = rast(:,2)>=extI_step_ends(i) & ...
                rast(:,2)<extI_step_ends(i+1);
            exc_fir_rate_first_half_vec(i) = ...
                sum(first_half_idx & local_fir_time_idx)* ...
                (t_unit/local_time_window_length)/(en/2);
            exc_fir_rate_second_half_vec(i) = ...
                sum(second_half_idx & local_fir_time_idx)* ...
                (t_unit/local_time_window_length)/(en/2);
            if in ~= 0
                inh_fir_rate_vec(i) = ...
                    sum(inh_idx & local_fir_time_idx)* ...
                    (t_unit/local_time_window_length)/((-1)*in);
            end
        end
    end

    % Plot
    plot_pop_firing_rate_vs_time = 0;
    plot_pop_firing_rate_vs_extI = 1;
    data_title_type = 4;

    if data_title_type == 4
        data_title = [num2str(en, 'e_n=%i, '), ...
            num2str(in, 'i_n=%i, ') ...
            num2str(SE, 'SE=%1.4f, ') num2str(SI, 'SI=%1.4f')];
    elseif data_title_type == 3
        data_title = [num2str(en, 'e_n=%1.0e, '), ...
            num2str(in, 'i_n=%1.0e, ') ...
            num2str(t_end, 't_e_n_d=%1.0e, ') num2str(dt, 'dt=%1.0e, ')...
            num2str(SE, 'SE=%1.4f, ')];
    elseif data_title_type == 2
        data_title = [num2str(en, 'e_n=%1.0e, ') ...
            num2str(t_end, 't_e_n_d=%1.0e, ') num2str(dt, 'dt=%1.0e, ')...
            num2str(SE, 'SE=%1.4f, ') num2str(Imin, 'extI=[%3i, ') ...
            num2str(Imax, '%3i] in ') num2str(nr_extI_steps, '%3i steps')];
    elseif data_title_type == 1
        data_title = [num2str(SE, 'SE=%1.4f, ') ...
            num2str(poisson_rate, 'Poisson: rate=%3.2f per sec, ') ...
            num2str(poisson_strength, 'strength=%3.2f')];
    end
        
    if plot_pop_firing_rate_vs_time
        figure,
        plot(extI_step_ends(2:end)/t_unit, exc_fir_rate_first_half_vec);
        hold on;
        plot(extI_step_ends(2:end)/t_unit, exc_fir_rate_second_half_vec);
        legend('Exc. 1/2 (w ext. input)', 'Exc. 2/2 (w/o ext. input)');
        if in ~= 0
            plot(extI_step_ends(2:end)/t_unit, inh_fir_rate_vec);
            legend('Exc. 1/2 (w ext. input)', ...
                'Exc. 2/2 (w/o ext. input)', 'Inh. (w/o ext. input)');
        end
        xlabel('t (in seconds)');
        ylabel('Firing rate (per sec)');
        title(data_title);
    end
    
    if plot_pop_firing_rate_vs_extI
        figure,
        plot(extI_step_values, exc_fir_rate_first_half_vec);
        hold on;
        plot(extI_step_values, exc_fir_rate_second_half_vec);
        legend('Exc. 1/2 (w ext. input)', 'Exc. 2/2 (w/o ext. input)');
        if in ~= 0
            plot(extI_step_values, inh_fir_rate_vec);
            legend('Exc. 1/2 (w ext. input)', ...
                'Exc. 2/2 (w/o ext. input)', 'Inh. (w/o ext. input)');
        end
        xlabel('ext_I');
        xlim([Imin Imax]);
        ylabel('Firing rate (per sec)');
        title(data_title);
    end
    
    % Save workspace
    if save_workspace
        filename = ['saved_workspaces/' 'LIF_test_' ...
            num2str(test_type_num) ...
            datestr(now,'_mmmdd_HHMMSS') '.mat'];
        % Exclude large variables
        save(filename, '-regexp', ['^(?!(' 'local_fir_time_idx|'...
            'rast|' 'first_half_idx|' 'second_half_idx|' 'inh_idx|' ...
            'extI_matrix', ...
            ')$).']);
    end
    
elseif test_type_num == 3
    %% HALF WITH EXTERNAL INPUTS, OTHER HALF WITHOUT
    % David's message:
    % Now for a network of, say, 400 hundred excitatory neurons of
    % all-to-all connections, suppose half of the neurons receive external
    % inputs, what is the gain curve for these two populations,
    % i.e. one half with external inputs, the other half without.
    
    % Input parameters
    t_unit = 1000;
    en = 400;
    in = 0;
    edge_prob = 1;
    t_end = 100;
    dt = 0.1;
    extI_vec = [8100:-500:7100]';
    nr_extI = length(extI_vec);
    rast{nr_extI} = [];
    
    % Model
    tSIM = tic;
    for i = 1:nr_extI
        [~, ~, ~, rast{i}] = LIF_network('ER', [en in], edge_prob, ...
            [t_end dt], 'external_spike_rate', 0, ...
            'display_input_output_summary', 0, 'plot_results', 0, ...
            'extI', [extI_vec(i)*ones(en/2,1); zeros(en/2,1)]);
        fprintf('[%3.2fm]: %i/%i completed.\n', toc(tSIM)/60, i, nr_extI);
    end
    
    % Crunch
    exc_fir_rate_first_half_vec = zeros(nr_extI, 1);
    exc_fir_rate_second_half_vec = zeros(nr_extI, 1);
    for i=1:nr_extI
        exc_fir_rate_first_half_vec(i) = ...
            sum(rast{i}(:,1)<en/2+1)*(t_unit/t_end)/(en/2);
        exc_fir_rate_second_half_vec(i) = ...
            sum(rast{i}(:,1)>en/2)*(t_unit/t_end)/(en/2);
    end
    figure,
    plot(extI_vec, exc_fir_rate_first_half_vec, 'x-');
    hold on;
    plot(extI_vec, exc_fir_rate_second_half_vec, 'x-');
    legend('Half with ext. input', 'Half w/o ext. input');
    xlabel('External input current (per sec)');
    ylabel('Firing rate (per sec)');
    
elseif test_type_num == 2
    %% CONSTANT INPUT TO A SINGLE NEURON
    % Input parameters
    t_unit = 1000;
    gL = 50/t_unit;
    en = 1;
    in = -1;
    edge_prob = 0;
    t_end = 10000;
    dt = 1;
    extI_vec = logspace(log10(50.00001), log10(100), 10)';
    % Model
    pfrm = run_LIF_network_with_multiple_values('ER', [en in], ...
        [t_end dt], 'network_params_vec', edge_prob, ...
        'input_rate_vec', 0, ...
        'extI_vec', extI_vec, 'plot_stuff', 0);
    % Crunch
    extI_theoretical = logspace(log10(min(extI_vec)), ...
        log10(max(extI_vec)), 1000);
    theoretical_firing_rate = -gL./log(1-gL./(extI_theoretical/t_unit))...
        *t_unit;
    
    sim_firing_rate = pfrm(:,1,1);
    figure,
    plot(extI_vec, sim_firing_rate, 'x');
    hold on;
    plot(extI_theoretical, theoretical_firing_rate);
    xlim([min(extI_vec)-1 max(extI_vec)+1]);
    legend('Simulation', 'Theory');
    title('Single Neuron - Constant input current');
    xlabel('External input current (per sec)');
    ylabel('Firing rate (per sec)');
    
elseif test_type_num == 1
    %% RK2 ERROR CONVERGENCE
    % Input parameters
    test_choice = 'fast';
    en = 2;
    in = -2;
    t_end = 100;   
    if strcmp(test_choice, 'ultra fast')
        d_reference = 1e-2;
        dt_vec = [0.1 0.2 0.4 0.8];
    elseif strcmp(test_choice, 'fast')
        d_reference = 1e-3;
        dt_vec = [0.01 0.02 0.04 0.08];
    elseif strcmp(test_choice, 'slow')
        d_reference = 1e-4;
        dt_vec = [0.00125 0.0025 0.005 0.01 0.02 0.04 0.08];
    end
    % Following are required for a deterministic simulation
    external_spike_type = 'Deterministic';
    external_spike_rate = 19999;
    type = 'ER';
    edge_prob = 1;
    SI = 0.05;
    track_stuff = 1;
    STDP = 'on';
    % Model
    n = en-in;
    nr_dt = length(dt_vec);
    % Evaluate the reference potential value
    v_ref =  LIF_network(type, [en in], edge_prob, ...
        [t_end d_reference], 'SI', SI, 'STDP', STDP, ...
        'external_spike_type', external_spike_type, ...
        'external_spike_rate', external_spike_rate, ...
        'track_sample_neurons', track_stuff, 'track_raster', track_stuff);
    v = zeros(n, nr_dt);
    for i = 1:nr_dt
        v(:,i) = LIF_network(type, [en in], edge_prob, ...
            [t_end dt_vec(i)], 'SI', SI, 'STDP', STDP, ...
            'external_spike_type', external_spike_type, ...
            'external_spike_rate', external_spike_rate, ...
        'track_sample_neurons', track_stuff, 'track_raster', track_stuff);
    end
    % Crunch
    addpath '../../../../Research/Sleep_Wake_Matlab_Codes_2013/test_rbm_data'
    
    % Take values of v(t) for the smallet dt as the true values, only
    % consider the first and the last entries; one excitatory and one
    % inhibitory
    true_vals = v_ref([1 end]);
    % Calculate errors
    errors = abs(v([1 end], :) - true_vals*ones(1,nr_dt));
    % Plot error vs dt and fit the regression line
    for i=1:2
        est_scaling_exponent_lin_reg(log(dt_vec), ...
            log(errors(i,:)), 'plot_results', 1);
    end
    
end

end