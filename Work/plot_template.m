
close all

figure
plot(t_test, y_test(1,:)); hold on;
plot(t_test, y_hat(1,:));

% change settings
plt = Plot(); % create a Plot object and grab the current figure

plt.XLabel = 'Time, t (s)';
plt.YLabel = 'Position, x (m)';
plt.XLim = [min(t_test) max(t_test)];
plt.YLim = [2 6.5];
plt.Colors = {                 % colors for different data sets
    [0.25,      0.25,       0.25]    
    [1,      0,       0]       
    };
plt.LineWidth = [2, 2];       % three line widths
plt.LineStyle = {'-', '--'}; % three line styles
plt.XMinorTick = 'off';
plt.YMinorTick = 'off';
plt.TickLength = [0.01, 0.01];
set(groot, 'defaultLegendInterpreter', 'latex'); % Set legend interpretter to Latex
plt.Legend = {'\textsf{Actual} $$x$$', '\textsf{Predicted} $$x$$'}; % legends
% plt.Legend = {'$$\theta$$', '$$\hat{\theta}$$'}; % legends

% Save? comment the following line if you do not want to save
sigma_str = strrep(num2str(sigma),'.','_');
train_str = num2str(N_train*Ts);
filename = strcat('HAVOK of x', ' train=', train_str, ' sig= ', sigma_str);
plt.export(strcat(filename,'.png'));  % Save in project folder
plt.export(strcat(filename,'.fig'));  % Save in project folder

% apply the settings
% setPlotProp(plt);