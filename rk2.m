function varargout = rk2(varargin)
%% Input parameters
dt = varargin{1};
v0 = varargin{2};
a0 = varargin{3};
a1 = varargin{4};
b0 = varargin{5};
b1 = varargin{6};
vI = varargin{7};
% ------------------------------------------------------------------------
%% Model
k1 = evaluate_k(v0, a0, b0);
k2 = evaluate_k(v0 + k1.*dt, a1, b1);
v1 = v0 + (k1+k2).*dt.*0.5;
%% Check outputs
if any(v1<vI)
    error('FO: Potential can never be less than vI! Consider smaller dt');
end
% ------------------------------------------------------------------------
%% Outputs
varargout{1} = v1;
end


% Old code in evolve_potential.m
% k1(idx) = evaluate_k(v0(idx), a0, b0(idx));
% k2(idx) = evaluate_k(v0(idx) + k1(idx).*dt, a1(idx), b1(idx));
% v1(idx) = v0(idx) + (k1(idx)+k2(idx)).*dt.*0.5;