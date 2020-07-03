
close all

figure
plot(t_test, y_test(1,:)); hold on;
plot(t_test, y_hat(1,:));

% change settings
plt = Plot(); % create a Plot object and grab the current figure

plt.XLabel = 'Time, t (s)';
plt.YLabel = 'Position, x (m)';
plt.Title = 'EKF model prediction of \x';
plt.XLim = [min(t_test) max(t_test)];
plt.Colors = {                 % colors for different data sets
    [0.25,      0.25,       0.25]    
    [1,      0,       0]       
    };
plt.LineWidth = [2, 2];       % three line widths
plt.LineStyle = {'-', '--'}; % three line styles
set(groot, 'defaultLegendInterpreter', 'latex'); % Set legend interpretter to Latex
plt.Legend = {'\textsf{Actual} $$x$$', '\textsf{Predicted} $$\theta$$'}; % legends
% plt.Legend = {'$$\theta$$', '$$\hat{\theta}$$'}; % legends

% Save? comment the following line if you do not want to save
plt.export('EKF_cartpend_predict_theta.png');  %???Where does it save?

% apply the settings
% setPlotProp(plt);