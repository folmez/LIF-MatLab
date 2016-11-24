function varargout = shelley_and_tao_2001(varargin)
%% Input parameters
dt = varargin{1};
firing_times = varargin{2};
a0 = varargin{3};
a1 = varargin{4};
b0 = varargin{5};
b1 = varargin{6};
vI = varargin{7}(1);
vL = varargin{7}(2);
% ------------------------------------------------------------------------
%% Check inputs
if (any(firing_times<0) || any(firing_times>dt)) && ~isempty(a0)
    error('FO: No firing time should lie outside [t, t+dt]!');
end
% ------------------------------------------------------------------------
%% Model
v0 = (vL - 0.5.*firing_times.*(b0 + b1 - a1.*b0.*dt)) ./ ...
        (1 + 0.5.*firing_times.* ((-1).*a0 + (-1).*a1 + a1.*a0.* dt));
v1 = rk2(dt, v0, a0, a1, b0, b1, vI);

%% Check outputs
if any(v1<vI)
    error('FO: Potential can never be less than vI!');
end
% ------------------------------------------------------------------------
%% Outputs
varargout{1} = v1;
end

% Old code in evolve_potential.m
%     v0(idx) = (vL - 0.5.*firing_times.* ...
%         (b0(idx) + b1(idx) - a1(idx).*b0(idx).*dt)) ./ ...
%         (1 + 0.5.*firing_times.* ...
%         ((-1).*a0(idx) + (-1).*a1(idx) + a1(idx).*a0(idx).* dt));
%     k1(idx) = evaluate_k(v0(idx), a0(idx), b0(idx));
%     k2(idx) = evaluate_k(v0(idx) + k1(idx).*dt, a1(idx), b1(idx));
%     v1(idx) = v0(idx) + (k1(idx)+k2(idx)).*dt.*0.5;
