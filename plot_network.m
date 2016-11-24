function plot_network(varargin)
%% Input parameters
A = varargin{1};
en = varargin{2};
in = varargin{3};

n = en + (-1)*in;
%% Plot
% Generate n angles 
angles = linspace(0, 2*pi, n+1);
angles(end) = [];

% Generate x- and y-coordinates
x = cos(angles);
y = sin(angles);

% Plot excitatory neurons
plot(x(1:en), y(1:en), 'r.', 'MarkerSize', 15); hold on;
% Plot inhibitory neurons
plot(x(en+1:n), y(en+1:n), 'b.', 'MarkerSize', 15);
% Write neuron numbers
for i = 1:n
    text(x(i)*1.2, y(i)*1.2, num2str(i));
end
% Plot edges
for i = 1:n
    for j = i+1:n
        if A(i,j)==1
            if i<=en && j<=en
                plot([x(i) x(j)], [y(i) y(j)], 'r');
            elseif i>en && j>en
                plot([x(i) x(j)], [y(i) y(j)], 'b');
            else
                plot([x(i) x(j)], [y(i) y(j)], 'k');
            end
        end
    end
end
    
xlim([-1.3 1.3]);
ylim([-1.3 1.3]);
end