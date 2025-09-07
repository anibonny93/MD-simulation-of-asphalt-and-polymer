% Data: Average Interaction Parameter for SARA Fractions Over Oxidation Time
data = [0.0096, 0.1144, 0.1501, 0.1969;
        0.0096, 0.1275, 0.2012, 0.2117;
        0.0096, 0.1364, 0.2115, 0.2583;
        0.0096, 0.1498, 0.2375, 0.2583;
        0.0096, 0.1587, 0.2902, 0.2613];

% Labels for oxidation days
oxidation_days = {'Day 0','Day 1','Day 5','Day 10','Day 30'};

% Create figure
figure;
hold on;

% Line plot (X-axis is now oxidation days)
plot(data, '-o', 'LineWidth', 1.5, 'MarkerSize', 8);
%bar(data, 'grouped');

% Colormap for better visibility
colormap(jet);

% Set x-axis labels as oxidation days
set(gca, 'XTick', 1:5, 'XTickLabel', oxidation_days, 'FontSize', 12);
set(gcf, 'Units', 'inches', 'Position', [1, 1, 10, 8]);

% Labels and title
xlabel('Oxidation Time (Days)', 'FontWeight', 'bold', 'FontSize', 14);
ylabel('Average Interaction Parameter', 'FontWeight', 'bold', 'FontSize', 14);
title('Average Interaction Parameter Over Oxidation Time', 'FontSize', 14, 'FontWeight', 'bold');

% Add legend for SARA fractions
legend({'Saturates', 'Asphaltenes', 'Aromatics', 'Resins'}, 'Location', 'NorthWest', 'FontSize', 12);

% Grid and box for clarity
grid off;
box on;

% Rotate x-axis labels for better readability if needed
xtickangle(45);

hold off;

