function varargout = evolve_potential(varargin)
%% Input parameters
t = varargin{1};
dt = varargin{2};
get0 = varargin{3};
git0 = varargin{4};
get1 = varargin{5};
git1 = varargin{6};
v0 = varargin{7};
en = varargin{8}(1);
in = varargin{8}(2);
A = varargin{9};
eneff = varargin{10};
ineff = varargin{11};
eori = varargin{12};
lft = varargin{13};
track_id = varargin{14};
extI = varargin{15}; 
spike_timing_dependent_plasticity = varargin{16};
SEE_matrix = varargin{17};
track_synaptic_strength = varargin{18};
vE = varargin{19}(1);
vI = varargin{19}(2);
vL = varargin{19}(3);
gL = varargin{19}(4);
vThresh = varargin{19}(5);
SEI = varargin{19}(6);
SI = varargin{19}(7);
tau_r_AMPA = varargin{20}(1);
tau_d_AMPA = varargin{20}(2);
tau_r_GABA = varargin{20}(3);
tau_d_GABA = varargin{20}(4);
tref = varargin{21};            % Refractory period
rs = varargin{22};              % Refractory state of each neuron
synaptic_strength_lb = varargin{23}(1);
synaptic_strength_ub = varargin{23}(2);
% ------------------------------------------------------------------------
display_track_progress_single_neuron = 0;
disable_warnings = 0;
n = en + (-1)*in;

%% Initialize display
if display_track_progress_single_neuron && t==0
    fprintf('\tID\tTime\tEvent\tv(t)\tge-tot\tgi-tot\n');
    fprintf('\t%i', track_id);
    fprintf('\t%3.2f', 0);
    fprintf('\t');
    fprintf('\t%1.4f', v0(track_id));
    fprintf('\t%1.3f', get0(track_id));
    fprintf('\t%1.3f\n', git0(track_id));
end

%% Modified RK2 for regular LIF evolution
% Initialization
v1 = vL*ones(n,1);

% Calculate RK2 parameters
a0 = get_alpha(get0, git0, gL);
a1 = get_alpha(get1, git1, gL);
b0 = get_beta(get0, git0, gL, vE, vI, vL, extI);
b1 = get_beta(get1, git1, gL, vE, vI, vL, extI);

% 1: Evolve potentials of neurons that are not in refractory period
% Find index set of neurons that are not inrefractory period
idx = rs==0;
% Use modified RK2 to evolve their potentials
v1(idx) = rk2(dt, v0(idx), a0(idx), a1(idx), b0(idx), b1(idx), vI);

% 2: Evolve potentials of neurons that exit refractory period in this time
% step using Shelley and Tao, 2001
% Find index set of neurons that exit refractory period
idx = rs==-1;
% Transpose their spike times according to the refractory period
firing_times = lft(idx) + tref - t;
% Use Shelley and Tao 2001 to evolve their potentials accordingly
v1(idx) = shelley_and_tao_2001(dt, firing_times, a0(idx), a1(idx), ...
    b0(idx), b1(idx), [vI vL]);
% Reset their refractory state
rs(idx) = 0;

% 3: Update refractory state of neurons whose refractory states are 1, i.e.
% neurons which will not fire in this time step, to -1 if they exit
% refractory period in the next time step
idx = (rs==1) & (lft+tref<=t+2*dt);
rs(idx) = -1;

% Display potential at t+dt
if display_track_progress_single_neuron
    fprintf('\t%i', track_id);
    fprintf('\t%3.2f', t+dt);
    fprintf('\t');
    fprintf('\t%1.4f', v1(track_id));
    fprintf('\t%1.3f', get1(track_id));
    fprintf('\t%1.3f\n', git1(track_id));
end

%% Calibrate firing events using Shelley and Tao, 2001
% Initialization
n = length(get0);
rast_increment = [];
genJump = zeros(n,1);
henJump = zeros(n,1);
ginJump = zeros(n,1);
hinJump = zeros(n,1);
firing_round_counter = 0;
SE1E2_track_increment = [];
while 1
    firing_round_counter = firing_round_counter+1;
    
    % Break out if no membrane potential is above the threshold
    if all(v1(rs==0)<=vThresh)
        break;
    end
    
    % Display error is a neuron fires twice in one time-step
    if firing_round_counter>1
        error('FO: Choose smaller dt! A neuron fires twice.');
    end
    
    % Find index set of neurons which fires at this time timestep and which
    % are also not in refractory during this time step
    idx = find(v1>vThresh);
    
    % Find the firing times in [0,dt]
    firing_times = (vThresh-v0(idx))./(v1(idx)-v0(idx))*dt;
    
    % Find subset of idx which exit refractory period in [0, dt]
    idx1 = idx(firing_times+tref<=dt);
    ft1 = firing_times(firing_times+tref<=dt)+tref;
    % Evaluate the post-firing potential using Eq. (12) from Shelley and
    % Tao, 2001 for neurons which exit refractory period in this time step
    v1(idx1) = shelley_and_tao_2001(dt, ft1, ...
        a0(idx1), a1(idx1), b0(idx1), b1(idx1), [vI vL]);

    % Find subset of idx which exit refractory period in (t+dt, t+2t] and
    % set their refractory state to -1
    idx2 = idx(dt<firing_times+tref & firing_times+tref<=2*dt);
    rs(idx2) = -1;
    
    % Find subset of idx which exit refractory period in (t+2t, inf) and 
    % set their refractory state to 1 and set their potential to reset
    % potential
    idx3 = idx(2*dt<firing_times+tref);
    rs(idx3) = 1;
    v1(idx3) = vL;
            
    % Update the last firing times in [t, t+dt]
    lft_prev = lft;
    lft(idx) = t + firing_times;
            
    % Display calibrated potential if tracked neuron has fired
    if display_track_progress_single_neuron && ismember(track_id, idx)
        fprintf('\t');
        fprintf('\t');
        fprintf('\tFIRED');
        fprintf('\t%1.4f\n', v1(track_id));
    end
    
    % Find number of neurons that fired during this time-step
    nr_firings = length(idx);
    
    % Evaluate post synaptic impact of the firing events
    for i = 1:nr_firings
        % Identify the fired neuron and its firing time
        id = idx(i);
        tfire = firing_times(i);
        % Find the indices of its neighbors
        id_nbrs = find(A(:,id));
        nr_nbrs = length(id_nbrs);
        hNew = zeros(nr_nbrs, 1);
        if eori(id)==1
            % Excitatory neuron firing
            tau_r = tau_r_AMPA;
            tau_d = tau_d_AMPA;
            % Exc. neuron (#id) firing into its exc. neighbors (#id_nbrs)
            temp = id_nbrs<en+1;
            hNew(temp) = (SEE_matrix(id, id_nbrs(temp))'/eneff)/tau_r_AMPA;
            % Firing into inhibitory neurons
            hNew(id_nbrs>en) = (SEI/eneff)/tau_r_AMPA;
            [gJump, hJump] = analytic(dt-tfire, tau_d, tau_r, 0, hNew);
            genJump(id_nbrs) = genJump(id_nbrs) + gJump;
            henJump(id_nbrs) = henJump(id_nbrs) + hJump;
        elseif eori(id)==-1 
            % Inhibitory neuron firing
            tau_r = tau_r_GABA;
            tau_d = tau_d_GABA;
            % Firing into excitatory neurons
            hNew(id_nbrs<en+1) = (SI/ineff)/tau_r_GABA;
            % Firing into inhibitory neurons
            hNew(id_nbrs>en) = (SI/ineff)/tau_r_GABA;
            [gJump, hJump] = analytic(dt-tfire, tau_d, tau_r, 0, hNew);
            ginJump(id_nbrs) = ginJump(id_nbrs) + gJump;
            hinJump(id_nbrs) = hinJump(id_nbrs) + hJump;        
        end
        % Display action potentials received by the tracked neuron
        if display_track_progress_single_neuron && ...
                ismember(track_id, id_nbrs)
            if eori(id)==1
                fprintf('\tRECEIVED EXCITATORY FIRING\n');
            elseif eori(id)==-1
                fprintf('\tRECEIVED INHIBITORY FIRING\n');
            end
        end
    end
    
    % Update E-E synpatic strengths using STDP
    if ~isempty(intersect(idx, 1:en))
        [SEE_matrix, SE1E2_track_increment] = update_SEE_matrix(...
            spike_timing_dependent_plasticity, ...
            SEE_matrix, idx, en, A, lft_prev, lft, ...
            track_synaptic_strength, ...
            [synaptic_strength_lb synaptic_strength_ub]);
    end  
    
    % Record firing times in a raster column
    rast_increment = [idx lft(idx)];
end

%% Output parameters
varargout{1} = v1;
varargout{2} = lft;
varargout{3} = rs;
varargout{4} = rast_increment;
varargout{5} = genJump;
varargout{6} = henJump;
varargout{7} = ginJump;
varargout{8} = hinJump;
varargout{9} = SEE_matrix;
varargout{10} = SE1E2_track_increment;
end