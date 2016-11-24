%% 
function varargout = delta_synaptic_strength(varargin)
% This synaptic strength update function follows (Morrison 2008, Sec. 4.1)

if strcmp(varargin{1}, 'test')
    test_delta_synaptic_strength;
else
    % Input arguments
    SEE_vec = varargin{1};
    delta_t_vec = varargin{2};
    
    % In-script arguments
    tau_plus = 10;
    tau_minus = 10;
    
    % Model
    % Find positive and negative firing time differences
    pos_delta_t_idx = delta_t_vec>0;
    neg_delta_t_idx = delta_t_vec<0;
    
    % Intiliaze the synaptic strength difference vector
    delta_ss = zeros(size(SEE_vec));
    
    % Calculate synaptic strengt differences
    delta_ss(pos_delta_t_idx) = F_plus(SEE_vec(pos_delta_t_idx)) .* ...
        exp(-abs(delta_t_vec(pos_delta_t_idx))/tau_plus);
    delta_ss(neg_delta_t_idx) = -F_minus(SEE_vec(neg_delta_t_idx)) .* ...
        exp(-abs(delta_t_vec(neg_delta_t_idx))/tau_minus);
    
    % Output arguments
    varargout{1} = delta_ss;
end

end

%% Test
function test_delta_synaptic_strength
n = 1000;
min_t = -100;
max_t = 100;
SEE_vec = ones(n, 1);
delta_t_vec = linspace(min_t, max_t, n)';

delta_ss = delta_synaptic_strength(SEE_vec, delta_t_vec);

figure,
idx = 1:n/2;
plot(delta_t_vec(idx), 100*(delta_ss(idx)./SEE_vec(idx)), 'r.');
hold on;
idx = n/2+1:n;
plot(delta_t_vec(idx), 100*(delta_ss(idx)./SEE_vec(idx)), 'g.');
plot([min_t max_t], [0 0], 'k:');
plot([0 0], [-40 100], 'k:');
ylabel('Relative change in synaptic strength (%)');
xlabel('Spike timing (ms)');
legend('pre fires after post => depression', ...
    'pre fires before post => potentiation', 'Location', 'Best');
end

%% 
function S_vec_out = F_plus(S_vec_in), S_vec_out = 1.00*S_vec_in; end
function S_vec_out = F_minus(S_vec_in), S_vec_out = 0.40*S_vec_in; end