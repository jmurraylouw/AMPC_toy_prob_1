%% Plot for group meeting 10 August
% MAE vs N_train for LS, HAVOK and EKF
% All 3 on 1 graph for each sigma

%% Plot MAE_x to compare algorithms
close all
figure(1); hold on; % MAE_x

model_list = {'EKF', 'LS', 'HAVOK'};
sigma = 0.1;

for i = 1:3
    model_name = model_list{i}
    
    % load results for model
    sig_str = strrep(num2str(sigma),'.','_'); % Convert sigma value to string
    save_file = ['Data\', model_name, '_N_train_vs_error', '_sig=', sig_str, '.mat'];

    load(save_file)    

    semilogy(N_train_saved/100, MAE_saved(1,:)); % MAE_x

end
hold off;

%Change settings
plt = Plot(); % create a Plot object and grab the current figure

% plt.XLabel = 'Time, t (s)';
% plt.YLabel = 'x (m) / the';
plt.FontSize = 12;
% plt.XLim = [t(N-N_test-N_train),  100];
plt.YLim = [0, 10];
plt.YScale = 'log';
% plt.Colors = {                 % colors for different data sets
%     [ 0.16,     0.44,    1.00 ] % 1 - u
%     [ 1.00,     0.50,    0.10 ] % 2 - x_train
%     [ 0.44,     0.00,    0.99 ] % 3 - theta_train
%     };
% [0.25, 0.25, 0.25]
% [ 0.16,     0.44,    1.00 ] blue 
% [ 0.00,     0.57,    0.00 ] green 
% [ 0.17,     0.17,    0.17 ] black
% [ 0.93,     0.00,    0.00 ] red
% [ 1.00,     0.50,    0.10 ] orange
% [ 0.44,     0.00,    0.99 ] purple
% [ 0.75,     0.00,    0.75 ] pink
% [ 0.50,     0.50,    0.50 ] gray

%                u    train     test      hat       split
%                1    2    3    4    5    6    7    8  
plt.LineWidth   = [1, 1, 1];
plt.LineStyle   = {'-', '-', '-'}; 
plt.XMinorTick  = 'off';
plt.YMinorTick  = 'off';
plt.TickLength  = [0.01, 0.01];
% set(groot, 'defaultLegendInterpreter', 'latex'); % Set legend interpretter to Latex
%                1    2    3    4    5    6    7    8  

plt.Legend = model_list;
% % Save? comment the following line if you do not want to save
sigma_str = strrep(num2str(sigma),'.','_');
filename = strcat('compare algos N_train vs error x', ' sig= ', sigma_str);
plt.export(strcat(filename,'.png'));  % Save in project folder
plt.export(strcat(filename,'.fig'));  % Save in project folder

%% Plot error for each algo at different sigma
N_train = 2000;
sigma_list = [0.1; 0.01;   0];
havok_x = [0.0333; 0.009;  0.0004357];
ekf_x   = [0.0433; 0.0421; 0.0237];
ls_x    = [0.0487; 0.0098; 0.0091];

figure(2); hold on;
plot(sigma_list, [havok_x, ekf_x, ls_x])

%Change settings
plt = Plot(); % create a Plot object and grab the current figure

% plt.XLabel = 'Time, t (s)';
% plt.YLabel = 'x (m) / the';
plt.FontSize = 12;
% plt.XLim = [t(N-N_test-N_train),  100];
% plt.YLim = [0, 10];
plt.YScale = 'log';
% plt.Colors = {                 % colors for different data sets
%     [ 0.16,     0.44,    1.00 ] % 1 - u
%     [ 1.00,     0.50,    0.10 ] % 2 - x_train
%     [ 0.44,     0.00,    0.99 ] % 3 - theta_train
%     };
% [0.25, 0.25, 0.25]
% [ 0.16,     0.44,    1.00 ] blue 
% [ 0.00,     0.57,    0.00 ] green 
% [ 0.17,     0.17,    0.17 ] black
% [ 0.93,     0.00,    0.00 ] red
% [ 1.00,     0.50,    0.10 ] orange
% [ 0.44,     0.00,    0.99 ] purple
% [ 0.75,     0.00,    0.75 ] pink
% [ 0.50,     0.50,    0.50 ] gray

%                u    train     test      hat       split
%                1    2    3    4    5    6    7    8  
plt.LineWidth   = [1, 1, 1];
plt.LineStyle   = {'.', '.', '.'}; 
plt.XMinorTick  = 'off';
plt.YMinorTick  = 'off';
plt.TickLength  = [0.01, 0.01];
% set(groot, 'defaultLegendInterpreter', 'latex'); % Set legend interpretter to Latex
%                1    2    3    4    5    6    7    8  

plt.Legend = model_list;
% % Save? comment the following line if you do not want to save
filename = strcat('compare algos sigma vs error (x)', ' N_train= ', num2str(N_train));
plt.export(strcat(filename,'.png'));  % Save in project folder
plt.export(strcat(filename,'.fig'));  % Save in project folder

%% Plot x only
figure(1);

plot(t_train, y_train(1,:)); % 1
hold on;
plot(t_test, y_test(1,:)); % 2
plot(t_test, y_hat(1,:)); % 3
plot(t(N-N_test).*[1,1], ylim); % 4    
hold off;

%Change settings
plt = Plot(); % create a Plot object and grab the current figure

% 1 - x_train
% 2 - x_test
% 3 - x_hat
% 4 - test/train split

plt.FontSize = 12;
plt.XLim = [t(N-N_test-N_train),  100];
% plt.YLim = [-2, 6];
plt.Colors = {                 % colors for different data sets
    [ 1.00,     0.50,    0.10 ] % 2 - x_train
    [ 0.17,     0.17,    0.17 ] % 4 - x_test
    [ 0.93,     0.00,    0.00 ] % 6 - x_hat
    [ 0.50,     0.50,    0.50 ] % 8 - test/train split
    };
% [0.25, 0.25, 0.25]
% [ 0.16,     0.44,    1.00 ] blue 
% [ 0.00,     0.57,    0.00 ] green 
% [ 0.17,     0.17,    0.17 ] black
% [ 0.93,     0.00,    0.00 ] red
% [ 1.00,     0.50,    0.10 ] orange
% [ 0.44,     0.00,    0.99 ] purple
% [ 0.75,     0.00,    0.75 ] pink
% [ 0.50,     0.50,    0.50 ] gray

%                train test  hat   split
%                2     4     6     8  
plt.LineWidth = [0.8,  1.2,  1.7,  1.5];
plt.LineStyle = {'-',  '-',  ':',  '-' };
plt.XMinorTick = 'off';
plt.YMinorTick = 'off';
plt.TickLength = [0.01, 0.01];
% set(groot, 'defaultLegendInterpreter', 'latex'); % Set legend interpretter to Latex

%             1    2    3    4    
plt.Legend = {'x_train', 'x_test', 'x_hat'}; 

% % Save? comment the following line if you do not want to save
sigma_str = strrep(num2str(sigma),'.','_');
train_str = num2str(N_train*Ts);
filename = strcat('HAVOK of x', ' train=', train_str, ' sig= ', sigma_str);
plt.export(strcat(filename,'.png'));  % Save in project folder
plt.export(strcat(filename,'.fig'));  % Save in project folder

%% Plot theta only
figure(1);


plot(t_train, y_train(2,:)); % 3
hold on;
plot(t_test, y_test(2,:)); % 5
plot(t_test, y_hat(2,:)); % 7
plot(t(N-N_test).*[1,1], ylim); % 8 
hold off;

%Change settings
plt = Plot(); % create a Plot object and grab the current figure

% 3 - theta_train
% 5 - theta_test
% 7 - theta_hat
% 8 - test/train split

% plt.XLabel = 'Time, t (s)';
% plt.YLabel = 'x (m) / the';
plt.FontSize = 12;
plt.XLim = [t(N-N_test-N_train),  100];
plt.Colors = {                 % colors for different data sets
    [ 0.44,     0.00,    0.99 ] % 3 - theta_train
    [ 0.00,     0.57,    0.00 ] % 5 - theta_test
    [ 0.16,     0.44,    1.00 ] % 7 - theta_hat
    [ 0.50,     0.50,    0.50 ] % 8 - test/train split
    };
% [0.25, 0.25, 0.25]
% [ 0.16,     0.44,    1.00 ] blue 
% [ 0.00,     0.57,    0.00 ] green 
% [ 0.17,     0.17,    0.17 ] black
% [ 0.93,     0.00,    0.00 ] red
% [ 1.00,     0.50,    0.10 ] orange
% [ 0.44,     0.00,    0.99 ] purple
% [ 0.75,     0.00,    0.75 ] pink
% [ 0.50,     0.50,    0.50 ] gray

%                train test hat  split
%                3     5    7    8  
plt.LineWidth = [0.4,  1.2, 1.7, 1.5];
plt.LineStyle = {'-',  '-', ':', '-' }; 
plt.XMinorTick = 'off';
% plt.YMinorTick = 'off';
plt.TickLength = [0.01, 0.01];
set(gca,'Layer','top'); % Set tick marks above plot layer

%                1    2    3    4    5    6    7    8  
plt.Legend = {':', '-', '-', '-'};

%% Save? comment the following line if you do not want to save
sigma_str = strrep(num2str(sigma),'.','_');
train_str = num2str(N_train*Ts);
filename = strcat('HAVOK of theta', ' train=', train_str, ' sig= ', sigma_str);
plt.export(strcat(filename,'.png'));  % Save in project folder
plt.export(strcat(filename,'.fig'));  % Save in project folder

