% Implentation of Hankel Alternative View Of Koopman
% Grid search of parameters
% Saves all the results for different parameter combinations

close all;
clear all;

% Search space
q_min = 30; % Min value of q in grid search
q_max = 300; % Max value of q in grid search
q_increment = 1; % Increment value of q in grid search

q_search = q_min:q_increment:q_max; % List of q parameters to search in

% Extract data
simulation_data_file = 'floating_pend_2D_data_3';
load(['Data/', simulation_data_file, '.mat']) % Load simulation data

u_data  = out.u.Data';
x_data  = out.x.Data';
measured_states = [1,2,3];
y_data  = x_data(measured_states,:); % Measurement data (x, z, theta)
t       = out.tout'; % Time

% Adjust for constant disturbance / mean control values
u_bar = mean(u_data,2); % Input needed to keep at a fixed point
% u_bar = [0; 6*9.81];
u_data  = u_data - u_bar; % Adjust for unmeasured input

% Testing data - Last 50 s is for testing and one sample overlaps training 
N_test = 2000; % Num of data samples for testing
x_test = x_data(:,end-N_test+1:end);
y_test = y_data(:,end-N_test+1:end); % One sample of testing data overlaps for initial condition
u_test = u_data(:,end-N_test+1:end);
t_test = t(:,end-N_test+1:end);

% Data dimentions
n = size(x_data,1); % number of states
m = size(y_data,1); % number of measurements
l = size(u_data,1); % number of inputs
Ts = t(2)-t(1);     % Sample time of data
N  = length(t);     % Number of data samples

% Add noise
rng('default');
rng(1); % Repeatable random numbers
sigma = 0.01; % Noise standard deviation
y_data_noise = y_data + sigma*randn(size(y_data));

% Training data - Last sample of training is first sample of testing
N_train = 4000; % Number of sampels in training data
y_train = y_data_noise(:,end-N_test-N_train+2:end-N_test+1); % Use noisy data
u_train = u_data(:,end-N_test-N_train+2:end-N_test+1);
t_train = t(:,end-N_test-N_train+2:end-N_test+1);
    
% Create empty results table
VariableTypes = {'int16', 'double'}; % q, MAE
VariableNames = {'q',     'MAE_mean'};
for i = 1:m % Mae column for each measured state
    VariableNames = [VariableNames, strcat('MAE_', num2str(i))];
    VariableTypes = [VariableTypes, 'double'];
end
Size = [length(q_search), length(VariableTypes)];

% Read previous results
sig_str = strrep(num2str(sigma),'.','_'); % Convert sigma value to string
results_file = ['Data/dmdc_results_', simulation_data_file, '_sig=', sig_str, '.mat'];

try
    load(results_file);
    results(~results.q,:) = []; % remove empty rows
    results = [results; table('Size',Size,'VariableTypes',VariableTypes,'VariableNames',VariableNames)];
    
catch
    disp('No saved results file')  
    
    results = table('Size',Size,'VariableTypes',VariableTypes,'VariableNames',VariableNames);
    emptry_row = 1; % Keep track of next empty row to insert results 
end

tic;

% Grid search
for q = q_search
    q

    if ~isempty(find(results.q == q, 1)) 
        continue % continue to next p if this combo has been searched before
    end

    w = N_train - q + 1; % num columns of Hankel matrix
    D = (q-1)*Ts; % Delay duration (Dynamics in delay embedding)

    % Create Hankel matrix with measurements
    Y = zeros(q*m,w); % Augmented state with delay coordinates [Y(k); Y(k-1*tau); Y(k-2*tau); ...]
    for row = 0:q-1 % Add delay coordinates
        Y(row*m+1:(row+1)*m, :) = y_train(:, row + (0:w-1) + 1);
    end

    % DMD of Y
    Y2 = Y(:, 2:end  );
    Y1 = Y(:, 1:end-1);

    YU = [Y1; u_train(:, q:end-1);]; % Combined matrix of Y and U, above and below
    AB = Y2*pinv(YU); % combined A and B matrix, side by side
    A  = AB(:,1:q*m); % Extract A matrix
    B  = AB(:,(q*m+1):end);

    A = stabilise(A,10);

    % Compare to testing data
    % Initial condition (last entries of training data)
    y_hat_0 = zeros(q*m,1);
    for row = 0:q-1 % First column of spaced Hankel matrix
        y_hat_0(row*m+1:(row+1)*m, 1) = y_train(:, end - ((q-1)+1) + row + 1);
    end

    % Run model
    Y_hat = zeros(length(y_hat_0),N_test); % Empty estimated Y
    Y_hat(:,1) = y_hat_0; % Initial condition
    for k = 1:N_test-1
        Y_hat(:,k+1) = A*Y_hat(:,k) + B*u_test(:,k);
    end

    y_hat = Y_hat(end-m+1:end, :); % Extract only non-delay time series (last m rows)

    % Vector of Mean Absolute Error on testing data
    MAE = sum(abs(y_hat - y_test), 2)./N_test; % For each measured state

    % Save results
    results(emptry_row,:) = [{q, mean(MAE)}, num2cell(MAE')]; % add to table of results
    emptry_row = emptry_row + 1; 

    save(results_file, 'results', 'emptry_row')

end % q

toc;

% Save results
results(~results.q,:) = []; % remove empty rows
save(results_file, 'results', 'emptry_row')

best_row = find(results.MAE_mean == min(results.MAE_mean));
best_results = results(best_row,:)

function A = stabilise(A_unstable,max_iterations)
    % If some eigenvalues are unstable due to machine tolerance,
    % Scale them to be stable
    A = A_unstable;
    count = 0;
    while (sum(abs(eig(A)) > 1) ~= 0)       
        [Ve,De] = eig(A);
        unstable = abs(De)>1; % indexes of unstable eigenvalues
        De(unstable) = De(unstable)./abs(De(unstable)) - 10^(-14 + count*2); % Normalize all unstable eigenvalues (set abs(eig) = 1)
        A = Ve*De/(Ve); % New A with margininally stable eigenvalues
        A = real(A);
        count = count+1;
        if(count > max_iterations)
            break
        end
    end
end
