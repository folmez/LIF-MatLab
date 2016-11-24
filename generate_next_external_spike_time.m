function varargout = generate_next_external_spike_time(varargin)
%%
generate_deterministic_times_no_matter_what = 0;
% -------------------------------------------------------------------------
%%
if strcmp(varargin{1}, 'test')
    test_generate_next_external_spike_time;
elseif strcmp(varargin{1}, 'Poisson') && ...
        ~generate_deterministic_times_no_matter_what
    last_spike_time = varargin{2};      % in ms
    poisson_spike_rate = varargin{3};   % in #/ms
    % ---------------------------------------------------------------------
    nr_requested_times = length(last_spike_time);
    next_spike_time = last_spike_time-...
        log(rand(nr_requested_times,1))/poisson_spike_rate;
    % ---------------------------------------------------------------------
    varargout{1} = next_spike_time;
elseif strcmp(varargin{1}, 'Deterministic') || ...
        generate_deterministic_times_no_matter_what
    last_spike_time = varargin{2};              % in ms
    deterministic_spike_rate = varargin{3};     % in #/ms
    % ---------------------------------------------------------------------
    next_spike_time = last_spike_time+1/deterministic_spike_rate;
    % ---------------------------------------------------------------------
    varargout{1} = next_spike_time;
end
end

function test_generate_next_external_spike_time
%%
n = 1e6;
nr_bins = 50;
rate = 1;
spike_times = generate_next_external_spike_time(zeros(n,1), rate);    
figure,
h = histogram(spike_times, nr_bins, 'Normalization', 'probability');
title(['Histogram of ' num2str(n) ...
    ' interarrival times with rate = ' num2str(rate)]);
xlabel('Spike time');
ylabel('PDF');
hold on;
% Plot theoretical expectation
bin_edges = h.BinEdges;
bin_centers = 0.5*(bin_edges(1:end-1)+bin_edges(2:end));
bin_left_ends = [0 bin_edges(2:end-1)];
bin_right_ends = [bin_edges(2:end-1) inf];
expected_bin_pdfs = -exp((-1)*rate*bin_right_ends)+...
    exp((-1)*rate*bin_left_ends);
plot(bin_centers, expected_bin_pdfs, 'r.-');
end