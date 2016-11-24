function varargout = evolve_ge_poisson_input(varargin)
%% Input parameters
t = varargin{1};
dt = varargin{2};
gep0 = varargin{3};
hep = varargin{4};
nr_rps = varargin{5};
lpst = varargin{6};
tau_r_AMPA = varargin{7};
tau_d_AMPA = varargin{8};
external_spike_type = varargin{9};
fS = varargin{10};
fR = varargin{11};
rate = varargin{12};
track_id = varargin{13};
% ------------------------------------------------------------------------
display_track_progress_single_neuron = 0;
use_daiwei_jump = 1;
% ------------------------------------------------------------------------
%% Check parameters
% lpst: Last Poisson Spike Time. This must be greater than or equal to for
% all neurons before the update begins
% FO (10-12-16): 50*eps is added to prevent machine epsilon size errors to
% spoil this comparison.
if any(lpst+50*eps<t)
    error(['FO: Last Poisson spike times must not be smaller ' ...
        'than the smaller endpoint of the timestep!']);
end
% ------------------------------------------------------------------------
%% Initialize
n = length(gep0);       % Network size
Tprev = ones(n,1)*t;    % Previous update time for each neuron
gep1 = gep0;            % Updated conductance for each neuron
% Initialize display
if display_track_progress_single_neuron && t==0
    fprintf('\tID\tTime\tEvent\t#spikes\tgep\thep\n');
end
% ------------------------------------------------------------------------
%% Evolve
while 1
    % Break out if there is no update
    if all(lpst>=t+dt), break; end
    % Find index set of the neurons that received a last Poisson spike time
    % in the current time step
    idx = lpst<t+dt;
    % Evolve g and h to the next spiking time
    [gep1(idx), hep(idx)] = analytic(lpst(idx)-Tprev(idx), ...
        tau_d_AMPA, tau_r_AMPA, gep1(idx), hep(idx));
    % Make the Poisson jump for the next spiking time
    if ~use_daiwei_jump
        hep(idx) = hep(idx)+fS/tau_d_AMPA;
    else
        hep(idx) = hep(idx)+fS/tau_r_AMPA;
    end
    % Increase number of received Poisson spikes by 1
    nr_rps(idx) = nr_rps(idx)+1;
    % Update previous time
    Tprev(idx) = lpst(idx);
    % Generate new spike times
    lpst(idx) = generate_next_external_spike_time(external_spike_type, ...
        lpst(idx), rate);
    % Display evolution if tracked neuron received a spike
    if display_track_progress_single_neuron && idx(track_id)
        fprintf('\t%i', track_id);
        fprintf('\t%3.2f', Tprev(track_id));
        fprintf('\tSpike');
        fprintf('\t%i', nr_rps(track_id));
        fprintf('\t%1.1e', gep1(track_id));
        fprintf('\t%1.1e\n', hep(track_id));
    end
end
% Last evolution from previous time (Tprev) to t+dt
[gep1, hep] = analytic((t+dt)-Tprev, ...
    tau_d_AMPA, tau_r_AMPA, gep1, hep);
% Display final evolution
if display_track_progress_single_neuron
    fprintf('\t%i', track_id);
    fprintf('\t%3.2f', t+dt);
    fprintf('\t');
    fprintf('\t');
    fprintf('\t%1.1e', gep1(track_id));
    fprintf('\t%1.1e\n', hep(track_id));
end
% ------------------------------------------------------------------------
%% Outputs
varargout{1} = gep1;
varargout{2} = hep;
varargout{3} = nr_rps;
varargout{4} = lpst;
end