function varargout = get_external_input_current(varargin)
%% Input parameters
extI_type = varargin{1};
t_unit = 1000;
plot_stuff = 0;

%% Model
step_end_times = 0;
step_fun_all_values = 0;
switch extI_type
    case 'time-dependent #1'
        % This time-dependent external input current is a step function:
        % gradual increase from 0 to 100 and then ramping down very slowly
        % to 0. (David: "... need to go up and down really slowly")
        
        % Input parameters - continued
        en = varargin{2}(1);
        in = varargin{2}(2);
        n = en + (-1)*in;
        t_end = varargin{3};
        dt = varargin{4};
        Imin = varargin{5}(1);      % Minimum external current
        Imax = varargin{5}(2);      % Maximum external current
        nr_steps = varargin{5}(3);  % # of steps in the functions (even)
        
        % Number of time steps of the LIF algorithm
        nr_time_steps = t_end/dt;
        % Step function's height increase at each jump
        step_height = (Imax-Imin)/(nr_steps-1);
        % Values of the step function on the left half
        left_half = Imin + step_height * ...
            round(linspace(0, (nr_steps-1), nr_time_steps/2));
        % Step function values
        step_fun_all_values = [left_half left_half(end:-1:1)];
        % External input matrix
        extI_matrix = ones(n,1)*step_fun_all_values;
        % Half the excitatory neurons should receive no external input
        extI_matrix(en/2+1:n, :) = 0;
        % End times of steps
        step_end_idx = sort(find(diff(extI_matrix(1,:)))+1);
        step_end_times = [0, step_end_idx*dt, nr_time_steps*dt]';
        % Step function individual step values
        step_fun_ind_values = [Imin step_fun_all_values(step_end_idx)]';
        
        % Plot
        if plot_stuff
            figure,
            plot((1:nr_time_steps)*dt/t_unit, extI_matrix(1,:), '.');
            hold on;
            plot(step_end_times(1:end-1)/t_unit, ...
                [Imin extI_matrix(1,step_end_idx)], 'o');      
            stairs(step_end_times/t_unit, ...
                extI_matrix(1,[1 step_end_idx end]), 'k');
            title('Time-dependent external input function #1');
            xlabel('t (in seconds)');         
            ylim([Imin-5 Imax+5]);
            
            figure,
            stairs(step_end_times/t_unit, ...
                extI_matrix(1,[1 step_end_idx end]), 'k');
            title('Time-dependent external input function #1');
            xlabel('t (in seconds)');         
            ylim([Imin-5 Imax+5]);
        end
        
        % Divide by time conversion factor
        extI_matrix = extI_matrix/t_unit;
        
    case 'constant'
        % Input parameters - continued
        extI  = varargin{2};
        n = varargin{3};
       
        if numel(extI) == 1
            % Turn the single constant input into a matrix
            extI_vector = ones(n, 1)*extI;
        elseif numel(extI) == n
            % Reshape external input current vector
            extI_vector = reshape(extI, [n, 1]);
        else
            error(['FO: External input current must either be a ' ...
                ' constant or an nx1 vector']);
        end
        
        % Divide by time conversion factor
        extI_vector = extI_vector/t_unit;
        
    case 'constant-this-case-is-discontinued-on-11-22-2016'
        % Input parameters - continued
        extI  = varargin{2};
        n = varargin{3};
        nr_time_steps = varargin{4};
        
        if numel(extI) == 1
            % Turn the single constant input into a matrix
            extI_matrix = ones(n, nr_time_steps)*extI;
        elseif numel(extI) == n
            % Reshape external input current vector
            extI = reshape(extI, [n, 1]);
            % Turn the vector into a matrix
            extI_matrix = extI*ones(1, nr_time_steps);
        else
            error(['FO: External input current must either be a ' ...
                ' constant or an nx1 vector']);
        end
    otherwise
        error('FO: Wrong external input type');
end

%% Outputs
switch extI_type
    case 'time-dependent #1'
        varargout{1} = extI_matrix;
        varargout{2} = step_end_times;
        varargout{3} = step_fun_ind_values;
    case 'constant'
        varargout{1} = extI_vector;
end
end