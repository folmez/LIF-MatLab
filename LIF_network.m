function varargout = LIF_network(varargin)
%% Network parameters
type = varargin{1};     % Type of network: 'ER'
en = varargin{2}(1);    % # of excitatory neurons (positive integer)
in = varargin{2}(2);    % # of inhibitory neurons (negative integer)
sn = 1;                 % # of sampled neurons from each population
params = varargin{3};   % Other network parameters
t_end = varargin{4}(1); % in milliseconds
dt = varargin{4}(2);    % in milliseconds

nr_pops = 2;            % # of populations: excitatory and inhibitory

display_input_output_summary = 1;
plot_results = 1;   
percent = 1.00;         % Fraction of results from the end to show on plot

track_raster = 1;
track_sample_neurons = 0;
track_synaptic_strength = 0;

external_spike_type = 'Poisson';    % 'Poisson' or 'Deterministic'

external_spike_rate = 4000;         % #/s
fS = 0.005;     % Strength of Poisson input
SE = 0.05;      % Excitatory synaptic strength 
SI = 0.5;       % Inhibitory synaptic strength

tref = 0;           % Refractory period

% Time constants
tau_r_AMPA = 1;
tau_d_AMPA = 5;
tau_r_GABA = 1;
tau_d_GABA = 5;

% External input parameters
extI_type = 'constant';
extI = 0;               % External input current
                        % 1x1 or 1xn

% STDP parameters                        
spike_timing_dependent_plasticity = 'off';
synaptic_strength_lb = 0.00;
synaptic_strength_ub = 0.10;
plot_subset_weight_distributions_at_every_progress_update = 1;
subset_size = 1000;
% ------------------------------------------------------------------------
i = 5;
while i<=length(varargin),
    switch varargin{i},
        case 'extI_matrix'
            extI_matrix = varargin{i+1};
        case 'extI'
            extI = varargin{i+1};
        case 'STDP'
            spike_timing_dependent_plasticity = varargin{i+1}{1};
            if length(varargin{i+1}) > 1
                synaptic_strength_lb = varargin{i+1}{2}(1);
                synaptic_strength_ub = varargin{i+1}{2}(2);
            end
        case 'external_spike_type', external_spike_type = varargin{i+1};
        case 'external_spike_rate', external_spike_rate = varargin{i+1};
        case 'external_spike_strength', fS = varargin{i+1};
        case 'rate_strengths'
            external_spike_rate = varargin{i+1}(1);
            fS = varargin{i+1}(2);
            SE = varargin{i+1}(3);
            SI = varargin{i+1}(4);
        case 'tau',
            tau_r_AMPA = varargin{i+1}(1);
            tau_d_AMPA = varargin{i+1}(2);
            tau_r_GABA = varargin{i+1}(3);
            tau_d_GABA = varargin{i+1}(4);
        case 'SE',                  SE = varargin{i+1};
        case 'SI',                  SI = varargin{i+1};
        case 'tref',                tref = varargin{i+1};
        case 'track_raster',        track_raster = varargin{i+1};
        case 'track_sample_neurons'
            track_sample_neurons = varargin{i+1};
        case 'track_synaptic_strength'
            track_synaptic_strength = varargin{i+1};
        case 'display_input_output_summary'
            display_input_output_summary = varargin{i+1};
        case 'plot_results'
            plot_results = varargin{i+1}(1);
            if length(varargin{i+1})>=2
                percent = varargin{i+1}(2);
            end
        otherwise,
            display(varargin{i});
            error('Unexpected inputs!!!');
    end
    i = i+2;
end
% ------------------------------------------------------------------------
%% Generate network
if in>0
    warning('Number of inhibitory neurons changed to a negative number!');
    in = -in;
end
n = en + (-1)*in;       % Total network size
network_params = {type, {en, in}, {params}};    % type - size - parameters
[A, eneff, ineff, ei_nbrs] = ...
    generate_LIF_network_adjacency_matrix(network_params);

%% Model parameters
t_unit = 1000;      % time conversion factor from ms to s
nr_time_steps = t_end/dt;

fR = 1;             % Ratio of Poisson input

vE = 14/3;          % Excitatory synaptic reversal potential
vI = -2/3;          % Inhibitory synaptic reversal potential
vL = 0;             % Reversal (Reset) potential for leakage
gL = 50/t_unit;     %
vThresh = 1;        % Normalized membrane potential threshold

%% Check inputs
if ~isinteger(int8(nr_time_steps))
    error('FO: Total time must be an integer multiple of dt!');
end

%% 

% External input current vector or matrix representing the external input
% current function
if strcmp(extI_type, 'constant')
    ext_vector = get_external_input_current(extI_type, extI, n);
elseif strcmp(extI_type, 'time-dependent #1')
    % extI_matrix is an n-by-nr_time_steps matrix whose entries hold the
    % external input current to each neuron at any given time step.
    extI_matrix = get_external_input_current(...
        extI_type, extI, n, nr_time_steps);
end

% Synaptic strength matrix representing strength between exc. neurons
if ~exist('SEE_matrix', 'var')
    % SEE_matrix is an en-by-en matrix, i.e. its number of rows and columns
    % are equal to the number of excitatory neurons in the model. Its
    % entries hold the synaptic strength values for each E-E synapse. The
    % entry at i-th row and j-th column, that is SEE_matrix(i,j), holds the
    % synaptic strength between E_i (presynaptic) and E_j (postsynaptic).
    SEE_matrix = ~eye(en) * SE; % presynaptic: E, postsynaptic: E
    SEI = SE;                   % presynaptic: E, postsynaptic: I
end

%% Initialization
% Neuron types in the network: first en neurons are excitatory (1) and the
% last in neurons are inhibitory (-1)                                
eori = [ones(en,1) ; (-1)*ones((-1)*in,1)];

lpst = zeros(n,1);  % Last external Poisson spike time into each neuron
gep1 = zeros(n,1);  % Poisson conductance at t_j+dt
hep = zeros(n,1);   
nr_rps = zeros(n,1);% # of received Poisson spikes

gen1 = zeros(n,1);  % Excitatory network conductance at t_j+dt
hen = zeros(n,1);

gin1 = zeros(n,1);  % Inhibitory network conductance at t_j+dt
hin = zeros(n,1);

v1 = zeros(n,1);    % Membrane potential at t_j+dt
lft = zeros(n,1);   % Last firing time of each neuron (i.e. v>vThresh)

rs = zeros(n,1);    % Refractory state of each neuron
                    %  0: Not in refractory state
                    % -1: In refractory period but will exit in the next
                    %     time step
                    %  1: In refractory period and will exit after the next
                    %     time step
warm_up_time = 0;

%% Initilization for tracking
% Raster matrix: [neuron ID ; Firing time]
if track_raster
    rast = zeros(nr_time_steps, 2);
    rast_count = 0;% #/ms
end

% Conductance and potential tracking
if track_sample_neurons
    gep_track = zeros(nr_pops*sn, nr_time_steps+1);
    hep_track = zeros(nr_pops*sn, nr_time_steps+1);
    nr_rps_track = zeros(nr_pops*sn, nr_time_steps+1);
    gen_track = zeros(nr_pops*sn, nr_time_steps+1);
    hen_track = zeros(nr_pops*sn, nr_time_steps+1);
    gin_track = zeros(nr_pops*sn, nr_time_steps+1);
    hin_track = zeros(nr_pops*sn, nr_time_steps+1);
    v_track = zeros(nr_pops*sn, nr_time_steps+1);
    % Select a random index set to track from each population
    if en ~= 0 && in ~= 0
        idx_set_track = union(datasample(1:en, sn, 'Replace', false), ...
            datasample(en+1:en+(-1)*in, sn, 'Replace', false));
        nr_sample_neurons = 2;
    elseif en ~= 0 && in == 0
        idx_set_track = sort(datasample(1:en, sn, 'Replace', false));
        nr_sample_neurons = 1;
    elseif en == 0 && in ~= 0
        idx_set_track = sort(datasample(en+1:en+(-1)*in, sn, ...
            'Replace', false));
        nr_sample_neurons = 1;
    end
    % Select a single neuron to track further
    single_track_id = datasample(idx_set_track, 1);
else
    single_track_id = 0;
end

% Synaptic strength tracking
if track_synaptic_strength
    SE1E2_track = zeros(2*nr_time_steps, 4); % time-id-E1E2-E2E1
    SE1E2_track_count = 0;
end

%% Model
% Turn rate from #/s to #/ms
rate = external_spike_rate/t_unit;  
% Last Poisson spike time into each neuron has to be initialized with the
% first spiking time in the beginning
lpst = generate_next_external_spike_time(external_spike_type, lpst, rate);
% Initialize progress tracking of the loop
progress_unit = 0.01;
current_progress = 0;
% Keep track of the whole simulations
tSIM = tic;
% Keep track display renewal
tDisplay = tic;
% Initialize weight distribution figure
if strcmp(spike_timing_dependent_plasticity, 'on') && ...
        plot_subset_weight_distributions_at_every_progress_update
    h1 = figure;
    hold on;
end
for i = 1:nr_time_steps
    % Initialize for evolution from t to t+dt
    t = (i-1)*dt;               % time
	gep0 = gep1;                % Poisson conductance at t
    gen0 = gen1;                % excitatory network conductance at t
    gin0 = gin1;                % inhibitory network conductance at t
    v0 = v1;                    % membrane potential at t
    % external input current
    if strcmp(extI_type, 'constant')
        extI = ext_vector;
    else
        extI = extI_matrix(:, i);
    end
    
    % Evolve Poisson conductance input to t+dt
    [gep1, hep, nr_rps, lpst] = evolve_ge_poisson_input(t, dt, ...
        gep0, hep, nr_rps, lpst, tau_r_AMPA, tau_d_AMPA, ...
        external_spike_type, fS, fR, rate, single_track_id);
    
    % Evolve network condtuctances and h to next time step
    [gen1, hen] = analytic(dt, tau_d_AMPA, tau_r_AMPA, gen0, hen);
    [gin1, hin] = analytic(dt, tau_d_GABA, tau_r_GABA, gin0, hin);
    
    % Update total conductances
    get0 = gen0 + gep0;
    get1 = gen1 + gep1;
    git0 = gin0;
    git1 = gin1;
    
    % Use total conductances to solve LIF and evolve membrane potential to
    % the next time step t+dt. Also, determine neurons that have reached
    % the threshold and evaluate how much network conductances need to be
    % jumped.
    [v1, lft, rs, rast_inc, genJump, henJump, ginJump, hinJump, ...
        SEE_matrix, SE1E2_track_inc] = evolve_potential(t, dt, ...
        get0, git0, get1, git1, ...
        v0, [en in], A, eneff, ineff, eori, lft, single_track_id, ...
        extI, spike_timing_dependent_plasticity, SEE_matrix, ...
        track_synaptic_strength, ...
        [vE vI vL gL vThresh SEI SI], ...
        [tau_r_AMPA tau_d_AMPA tau_r_GABA tau_d_GABA], ...
        tref, rs, [synaptic_strength_lb synaptic_strength_ub]);
    % Evolve network conductances
    gen1 = gen1+genJump;
    hen = hen+henJump;
    gin1 = gin1+ginJump;
    hin = hin+hinJump;
    
    % Update the raster matrix and raster count
    if track_raster && ~isempty(rast_inc)
        new_rast_count = rast_count+size(rast_inc, 1);
        rast(rast_count+1:new_rast_count, :) = rast_inc;
        rast_count = new_rast_count;
    end
    
    % Record data for tracked neurons
    if track_sample_neurons
        gep_track(:, i+1) = gep1(idx_set_track);
        hep_track(:, i+1) = hep(idx_set_track);
        nr_rps_track(:, i+1) = nr_rps(idx_set_track);
        gen_track(:, i+1) = gen1(idx_set_track);
        hen_track(:, i+1) = hen(idx_set_track);
        gin_track(:, i+1) = gin1(idx_set_track);
        hin_track(:, i+1) = hin(idx_set_track);
        v_track(:, i+1) = v1(idx_set_track);
    end
    
    % Update synaptic strength tracker
    if track_synaptic_strength && ~isempty(SE1E2_track_inc)
        new_SE1E2_track_count = SE1E2_track_count+size(SE1E2_track_inc,1);
        SE1E2_track(SE1E2_track_count+1:new_SE1E2_track_count, :) = ...
            SE1E2_track_inc;
        SE1E2_track_count = new_SE1E2_track_count;
    end
    
    % Display simulation progress
    if toc(tDisplay) > 6
        % Renew display clock
        tDisplay = tic;
        fprintf('%3.2f%% is completed in %3.2f minutes...\n', ...
            100*i/nr_time_steps, toc(tSIM)/60);
        current_progress = current_progress+progress_unit;

        % Plot synaptic strength histogram
        if strcmp(spike_timing_dependent_plasticity, 'on') && ...
                plot_subset_weight_distributions_at_every_progress_update
            SEE_vector = SEE_matrix(logical(A(1:en, 1:en)));           
            if length(SEE_vector) > subset_size
                SEE_vector = SEE_vector(randperm(end, subset_size));
                fig_tit = ['Histogram of a subset of weights (#=' ...
                num2str(subset_size) ')'];
            else
                fig_tit = ['Histogram of a subset of weights (#=' ...
                num2str(length(SEE_vector)) ')'];
            end
            figure(h1),
            histogram(SEE_vector, 'Normalization', 'countdensity', ...
                'DisplayStyle', 'stairs');
            xlabel('Synaptic strength');
            ylabel('Count density');
            title(fig_tit);
            drawnow;
        end       
    end
end

%% Calculate population statistics
exc_poisson_input = sum(nr_rps(1:en))*t_unit/t_end/en;
inh_poisson_input = sum(nr_rps(en+1:n))*t_unit/t_end/(-in);
tot_poisson_input = sum(nr_rps)*t_unit/t_end/n;

% Raster
if track_raster
    % Delete unused raster entries
    rast(rast_count+1:end, :) = [];
    if warm_up_time>t_end
        warning('FO: Warm up time is longer than total time.');
    end
    idx = rast(:,2)>warm_up_time;
    exc_firing_rate = sum(rast(idx,1)<en+1)*...
        (t_unit/(t_end-warm_up_time))/en;
    inh_firing_rate = sum(rast(idx,1)>en)*...
        (t_unit/(t_end-warm_up_time))/(-in);
    tot_firing_rate = sum(idx)*...
    (t_unit/(t_end-warm_up_time))/n;
else
    exc_firing_rate = [];
    inh_firing_rate = [];
    tot_firing_rate = [];
end

% Synaptic strength
if track_synaptic_strength
    SE1E2_track(SE1E2_track_count+1:end, :) = [];
end

% Record synaptic strengths of the existing EE synapses
if strcmp(spike_timing_dependent_plasticity, 'on')
    SEE_vector = SEE_matrix(logical(A(1:en, 1:en)));
end

avg_exc_v = mean(v1(1:en));
avg_inh_v = mean(v1(en+1:n));
avg_v = mean(v1);

avg_exc_gep = mean(gep1(1:en));
avg_inh_gep = mean(gep1(en+1:n));
avg_gep = mean(gep1);

avg_exc_gen = mean(gen1(1:en));
avg_inh_gen = mean(gen1(en+1:n));
avg_gen = mean(gen1);

avg_exc_gin = mean(gin1(1:en));
avg_inh_gin = mean(gin1(en+1:n));
avg_gin = mean(gin1);

%% Display input and output summary
if display_input_output_summary
    fprintf('               ----------------------------------\n');
    fprintf('               | LIF network simulation summary |\n');
    fprintf('               ----------------------------------\n');
    % Input parameters
    fprintf('\t               Network: %s (%1.4f)\n', type, params);
    fprintf('\t          [#exc, #inh]: [%i, %i]\n', en, in);
    fprintf('\t                  Neff: [%1.2f, %1.2f]\n', eneff, ineff);
    fprintf('\t       Refractory time: %1.2f ms\n', tref);
    if strcmp(spike_timing_dependent_plasticity, 'off')
        fprintf('\t     Cortical strength: [%1.2f, %1.2f]\n', SE, SI);
    end
    fprintf('\t                  STDP: %s\n', ...
        spike_timing_dependent_plasticity);
    fprintf('\t                    dt: %1.2f ms\n', dt);
    fprintf('\t                  Time: %i ms\n', t_end);
    fprintf('\t               Warm-up: %i ms\n', warm_up_time);
    fprintf('\n');
    if length(unique(extI))==1
        fprintf('\t    Ext. input current: %1.2f per sec\n', ...
            extI(1)*t_unit);
    else
        fprintf('\t    Ext. input current: %ix1 vector', n);
        fprintf(' in [%1.2f, %1.2f] per sec\n', ...
            min(min(extI_matrix))*t_unit, max(max(extI_matrix))*t_unit);
    end
    fprintf('\t       Ext. input type: %s', external_spike_type);
    fprintf(' (rate: %4.2f per sec, ', external_spike_rate);
    fprintf(', strength: %1.3f)\n', fS);
    % Output parameters
    fprintf('\n');
    fprintf('\t Poisson input to exc.: %4.2f per sec\n', ...
        exc_poisson_input);
    fprintf('\t Poisson input to inh.: %4.2f per sec\n', ...
        inh_poisson_input);
    fprintf('\t   Total Poisson input: %4.2f per sec\n', ...
        tot_poisson_input);
    fprintf('\n');
    fprintf('\tExcitatory firing rate: %4.3f per sec\n', exc_firing_rate);
    fprintf('\tInhibitory firing rate: %4.3f per sec\n', inh_firing_rate);
    fprintf('\t     Total firing rate: %4.3f per sec\n', tot_firing_rate);
    fprintf('\n');
    fprintf('\tAverage exc. potential: %1.3f\n', avg_exc_v);
    fprintf('\tAverage inh. potential: %1.3f\n', avg_inh_v);
    fprintf('\t     Average potential: %1.3f\n', avg_v);
    fprintf('\n');
    fprintf('\t      Average exc. gep: %2.3f\n', avg_exc_gep);
    fprintf('\t      Average inh. gep: %2.3f\n', avg_inh_gep);
    fprintf('\t           Average gep: %2.3f\n', avg_gep);
    fprintf('\n');
    fprintf('\t      Average exc. gen: %2.3f\n', avg_exc_gen);
    fprintf('\t      Average inh. gen: %2.3f\n', avg_inh_gen);
    fprintf('\t           Average gen: %2.3f\n', avg_gen);
    fprintf('\n');
    fprintf('\t      Average exc. gin: %2.3f\n', avg_exc_gin);
    fprintf('\t      Average inh. gin: %2.3f\n', avg_inh_gin);
    fprintf('\t           Average gin: %2.3f\n', avg_gin);
end

%% Plot sample
if plot_results
    if track_sample_neurons
        plot_inhibitory_stats = 1;
        tt = t_end*(1-percent):dt:t_end;
        ii = 1+(nr_time_steps*(1-percent):nr_time_steps);      
        figure,
        for k = 1:nr_sample_neurons
            subplot(3,2,k);
            plot(tt, gep_track(k,ii)); hold on;
            plot(tt, gen_track(k,ii));
            if plot_inhibitory_stats, plot(tt, gin_track(k,ii)); end
            xlim([tt(1) tt(end)]);
            subplot(3,2,k+2);
            plot(tt, hep_track(k,ii)); hold on;
            plot(tt, hen_track(k,ii));
            if plot_inhibitory_stats, plot(tt, hin_track(k,ii)); end
            xlim([tt(1) tt(end)]);
            subplot(3,2,k+4);
            plot(tt, v_track(k,ii));
            xlabel('ms');
            legend(['v(t) for #' num2str(idx_set_track(k))], ...
                'Location', 'Best');
            xlim([tt(1) tt(end)]);
            ylim([vI vThresh])
            hold on;
            plot([tt(1) tt(end)], [vL vL], 'r--');
        end
        subplot(3,2,2);
        if plot_inhibitory_stats, legend('g_e_p', 'g_e_n', 'g_i_n');
        else legend('g_e_p', 'g_e_n');
        end
        subplot(3,2,4);
        if plot_inhibitory_stats, legend('h_e_p', 'h_e_n', 'h_i_n');
        else legend('h_e_p', 'h_e_n');
        end
        subplot(3,2,6);
    end
    
    % Plot raster
    plot_network_connections = 0;
    if track_raster
        rast_exc_ii = (rast(:,2)>=t_end*(1-percent)) & rast(:,1)<en+1;
        rast_inh_ii = (rast(:,2)>=t_end*(1-percent)) & rast(:,1)>en;
        figure,
        if plot_network_connections
            subplot(1,2,1);
        end
        plot(rast(rast_exc_ii,2), rast(rast_exc_ii,1), 'r.', ...
            'MarkerSize', 2);
        hold on;
        plot(rast(rast_inh_ii,2), rast(rast_inh_ii,1), 'b.', ...
            'MarkerSize', 2);
        ylim([0 n+1]);
        xlim([t_end*(1-percent) t_end]);
        title('Raster plot');
        xlabel('ms');
        ylabel('Neuron ID');
        set(gca,'YDir','reverse')
        box on;
        
        if plot_network_connections
            subplot(1,2,2);
            plot_network(A, en, in);
            %
            %         sub(gca,'YDir','reverse')
            %         xlaplot(1,2,2);
            %         plot(ei_nbrs(:,1), 1:n, 'r.');
            %         hold on;
            %         plot(ei_nbrs(:,2), 1:n, 'b.');
            %         ylim([0 n+1]);
            %         setbel('# of nbr');
        end
    end
    drawnow;
    
    % Plot synaptic strengths E1->E2 and E2->E1
    if track_synaptic_strength && ...
            strcmp(spike_timing_dependent_plasticity, 'on')
        figure,
        plot(SE1E2_track(:,1), SE1E2_track(:,3));
        hold on;
        plot(SE1E2_track(:,1), SE1E2_track(:,4));
        xlabel('ms');
        ylabel('Synaptic strength');
        legend('E1->E2', 'E2->E1', 'Location', 'Best');
    end
    
    % Plot synaptic strength histogram
    if strcmp(spike_timing_dependent_plasticity, 'on')
        figure,
        histogram(SEE_vector, 'normalization', 'probability');
        xlabel('Synaptic strength');
    end
end

%% Output arguments
varargout{1} = v1;
varargout{2} = A;

if track_raster
    i = 3;
    varargout{i}(1) = tot_firing_rate;
    varargout{i}(2) = exc_firing_rate;
    varargout{i}(3) = inh_firing_rate;
    varargout{i+1} = rast;
end

if strcmp(spike_timing_dependent_plasticity, 'on')
    i = 3;
    varargout{i} = SEE_vector;
end

end