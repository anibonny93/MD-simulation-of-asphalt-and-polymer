% Data: Average Interaction Parameter for SARA Fractions Over Oxidation Time
data = [0.0096, 0.114, 0.103, 0.198;
        0.0096, 0.127, 0.213, 0.234;
        0.0096, 0.136, 0.218, 0.235;
        0.0096, 0.157, 0.254, 0.250;
        0.0096, 0.160, 0.291, 0.248];

data_stacked = data';
% Labels for oxidation days
xlabels = {'Saturates','Asphaltenes','Aromatics','Resins'};
legendlabels = {'Day 0', 'Day 1', 'Day 5', 'Day 10', 'Day 30'};

figure;
hold on;

bar(data_stacked, 'hist');
colormap(jet);

% Set x-axis labels as oxidation days
set(gca, 'XTick', 1:4, 'XTickLabel', xlabels, 'FontSize', 12);
set(gcf, 'Units', 'inches', 'Position', [1, 1, 10, 8]);

% Labels and title
xlabel('SARA Fractions', 'FontWeight', 'bold', 'FontSize', 14);
ylabel('Average Interaction Parameter (\chi)', 'FontWeight', 'bold', 'FontSize', 14);
title('Average Interaction Parameter Over Oxidation Time', 'FontSize', 14, 'FontWeight', 'bold');

% Add legend for SARA fractions
legend(legendlabels, 'Location','NorthWest');
legend({'Saturates', 'Asphaltenes', 'Aromatics', 'Resins'}, 'Location', 'NorthWest', 'FontSize', 12);

% Grid and box for clarity
grid on;
box on;

% Rotate x-axis labels for better readability if needed
xtickangle(45);

hold off;

