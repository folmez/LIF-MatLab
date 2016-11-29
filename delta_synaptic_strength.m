function varargout = delta_synaptic_strength(varargin)
% This synaptic strength update function follows (Morrison 2008, Sec. 4.1)

if strcmp(varargin{1}, 'test')
    which_type = varargin{2};
    test_delta_synaptic_strength(which_type);
else
    % Input arguments
    SEE_vec = varargin{1};
    delta_t_vec = varargin{2};
    % 'Bi1998':
    % 'Gutig2003':  See second page (pg. 3698) of Gutig et. al. 2003
    if length(varargin)>2
        which_type = varargin{3};
    else
        which_type = 'Gutig2003';
    end
    
    % Find positive and negative firing time differences
    pos_idx = delta_t_vec>0;
    neg_idx = delta_t_vec<0;
    
    % Intiliaze the synaptic strength difference vector
    delta_ss = zeros(size(SEE_vec));
    
    % Calculate synaptic strengt differences
    delta_ss(pos_idx) = F(which_type, 'plus', SEE_vec(pos_idx)) ...
        .* K(which_type, 'plus', delta_t_vec(pos_idx));
    delta_ss(neg_idx) = (-1)*F(which_type, 'minus', SEE_vec(neg_idx)) ...
        .* K(which_type, 'minus', delta_t_vec(neg_idx));
    
    % Check output
    if ~isreal(delta_ss)
        error('FO: Synaptic strength change is complex!');
    end
    
    % Output arguments
    varargout{1} = delta_ss;
end

end

%% Test
function test_delta_synaptic_strength(which_type)
n = 1000;
min_t = -100;
max_t = 100;
SEE_vec = ones(n, 1)*0.05;
delta_t_vec = linspace(min_t, max_t, n)';

delta_ss = delta_synaptic_strength(SEE_vec, delta_t_vec, which_type);

figure,
% Plot depression in red color on the left half
idx = 1:n/2;
plot(delta_t_vec(idx), 100*(delta_ss(idx)./SEE_vec(idx)), 'r.');
min_y_val_on_plot = min(100*(delta_ss(idx)./SEE_vec(idx)));
hold on;
% Plot potentiation in green color on the right half
idx = n/2+1:n;
plot(delta_t_vec(idx), 100*(delta_ss(idx)./SEE_vec(idx)), 'g.');
max_y_val_on_plot = max(100*(delta_ss(idx)./SEE_vec(idx)));
% Plot dotted lines along the axes
plot([min_t max_t], [0 0], 'k:');
plot([0 0], [min_y_val_on_plot max_y_val_on_plot], 'k:');

xlabel('Spike timing (ms)');
ylabel('Relative change in synaptic strength (%)');
legend('pre fires after post => depression', ...
    'pre fires before post => potentiation', 'Location', 'Best');
title(which_type);
end

%%
function S_vec_out = F(which_type, pm, S_vec_in)
switch which_type
    case 'Bi1998'
        if strcmp(pm, 'plus')
            S_vec_out = 1.00*S_vec_in;
        elseif strcmp(pm, 'minus')
            S_vec_out = 0.40*S_vec_in;
        end
    case 'Gutig2003'
        alpha = 1.05;  % Scale of asym. between potentiation and depression
        mu = 0.10;
        lambda = 0.05;  % Learning rate
        if strcmp(pm, 'plus')
            S_vec_out = (1-S_vec_in).^mu;
        elseif strcmp(pm, 'minus')
            S_vec_out = alpha*(S_vec_in.^mu);
        end
        S_vec_out = lambda*S_vec_out;
end
end

%%
function temporal_filter = K(which_type, pm, delta_t_vec)
switch which_type
    case 'Bi1998'
        tau_plus = 10;
        tau_minus = 10;
        if strcmp(pm, 'plus')
            temporal_filter = exp((-1)*abs(delta_t_vec)/tau_plus);
        elseif strcmp(pm, 'minus')
            temporal_filter = exp((-1)*abs(delta_t_vec)/tau_minus);
        end
    case 'Gutig2003'
        tau_plus = 20;
        tau_minus = 20;
        if strcmp(pm, 'plus')
            temporal_filter = exp((-1)*abs(delta_t_vec)/tau_plus);
        elseif strcmp(pm, 'minus')
            temporal_filter = exp((-1)*abs(delta_t_vec)/tau_minus);
        end
end
end