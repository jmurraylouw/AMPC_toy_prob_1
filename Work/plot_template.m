
%% Plot data vs model
figure(1);
plot(t, u_data); % 1
hold on;
plot(t_train, y_train); % 2,3 
plot(t_test, y_test); % 4,5
plot(t_test, y_hat); % 6,7 % Plot only non-delay coordinate
plot(t(N-N_test).*[1,1], ylim); % 8    
hold off;

%% Change settings
plt = Plot(); % create a Plot object and grab the current figure

% 1 - u
% 2 - x_train
% 3 - theta_train
% 4 - x_test
% 5 - theta_test
% 6 - x_hat
% 7 - theta_hat
% 8 - test/train split

% plt.XLabel = 'Time, t (s)';
% plt.YLabel = 'x (m) / the';
plt.FontSize = 15
plt.XLim = [t(N-N_test-N_train),  100];
plt.YLim = [-2, 6];
plt.Colors = {                 % colors for different data sets
    [ 0.16,     0.44,    1.00 ] % 1 - u
    [ 1.00,     0.50,    0.10 ] % 2 - x_train
    [ 0.44,     0.00,    0.99 ] % 3 - theta_train
    [ 0.17,     0.17,    0.17 ] % 4 - x_test
    [ 0.00,     0.57,    0.00 ] % 5 - theta_test
    [ 0.93,     0.00,    0.00 ] % 6 - x_hat
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

%                u    train     test      hat       split
%                1    2    3    4    5    6    7    8  
plt.LineWidth = [1.0, 0.8, 0.8, 1.2, 1.2, 1.7, 1.7, 1.5];       % three line widths
plt.LineStyle = {':', '-', '-', '-', '-', ':', ':', '-' }; % three line styles
plt.XMinorTick = 'off';
plt.YMinorTick = 'off';
plt.TickLength = [0.01, 0.01];
% set(groot, 'defaultLegendInterpreter', 'latex'); % Set legend interpretter to Latex
% plt.Legend = {'\textsf{Actual} $$x$$', '\textsf{Predicted} $$x$$'}; % legends
% plt.Legend = {'$$\theta$$', '$$\hat{\theta}$$'}; % legends

% Save? comment the following line if you do not want to save
sigma_str = strrep(num2str(sigma),'.','_');
train_str = num2str(N_train*Ts);
filename = strcat('HAVOK of x', ' train=', train_str, ' sig= ', sigma_str);
plt.export(strcat(filename,'.png'));  % Save in project folder
plt.export(strcat(filename,'.fig'));  % Save in project folder

