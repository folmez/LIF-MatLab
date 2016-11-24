function varargout = analytic(varargin)
%% Input parameters
dt = varargin{1};       % Timestep
tau_d = varargin{2};    % Decay constant
tau_r = varargin{3};    % Rise constant
g = varargin{4};        % Conductance (may be a column vector)
h = varargin{5};        % Conductance (may be a column vector)

%% Model
etd = exp((-1)*dt/tau_d);
etr = exp((-1)*dt/tau_r);
c = tau_r/(tau_d-tau_r).*(etd-etr);
g = g .* etd + c .* h;
h = h .* etr;

%% Outputs
varargout{1} = g;
varargout{2} = h;
end