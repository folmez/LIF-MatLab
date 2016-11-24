lambda = 4e5;
f = 0.005;
tau_r = 1;
tau_d = 5;
vE = 14/3;      % Excitatory synaptic reversal potential
vI = -2/3;      % Inhibitory synaptic reversal potential
vL = 0;         % Reversal (Reset) potential for leakage
gL = 50/1000; % 
vR = vL;

dt = 1000/lambda;
h0_left = (f/tau_r)/(exp(dt/tau_r)-1);
h0_right = h0_left+f/tau_r;

% Define the h function
h = @(t) h0_right*exp(-t/tau_r);

C = h0_right*tau_r/(tau_r-tau_d)*(exp(-dt/tau_r)-1)/(1-exp(-dt/tau_d));

% Define the g function (Poisson conductace)
g = @(t) C*exp(-t/tau_d) + tau_r/(tau_r-tau_d)*h0_right*exp(-t/tau_r);


g0 = g(0);
g1 = g(dt);

% Find v'(0) from the LIF equation
v_prime_at_0 = g0*vE + gL*vR;
% Find the time required for the membrane potential to reach 1
T1 = 1/v_prime_at_0;
% Find, approximately, v'(1)
T1_converted = (T1/dt-floor(T1/dt))*dt;
gep_at_T1 = g(T1_converted);
v_prime_at_T1 = -(gL+gep_at_T1)*1 + (gL*vR+gep_at_T1*vE);
% Find the time required for the membrane potential to reach 1
T2 = 1/v_prime_at_T1;
% Take the average of the two times
T = mean([T1 T2]);

T

fr= 1/T*1000;

fr



% 
% g0 = C + tau_r/(tau_r-tau_d)*h0_right;
% g1 = C*exp(-dt/tau_d) + tau_r/(tau_r-tau_d)*h0_right*exp(-dt/tau_r);
