% Data for AAA1 system
days_of_oxidation_AAA1 = [0, 1, 5, 10, 30];
oxygen_content_AAA1 = [0.39, 3.19, 5.04, 7.8, 9.4];
chi_value_AAA1 = [0.1804, 0.2085, 0.2269, 0.2376, 0.2567];

% Data for AAK1 system
days_of_oxidation_AAK1 = [0, 1, 5, 10, 30];
oxygen_content_AAK1 = [0.40, 3.29, 5.25, 8.06, 9.73];
chi_value_AAK1 = [0.1811, 0.2091, 0.2298, 0.2415, 0.2601];

% Data for AAM1 system
days_of_oxidation_AAM1 = [0, 1, 5, 10, 30];
oxygen_content_AAM1 = [0.41, 3.24, 4.29, 7.22, 8.15];
chi_value_AAM1 = [0.1248, 0.1564, 0.1883, 0.2010, 0.2257];

% Create the line plot for AAA1 system

figure;


%plot(days_of_oxidation_AAA1, chi_value_AAA1, 'bo-', 'MarkerFaceColor', 'b','LineWidth',3.5);
%plot(oxygen_content_AAA1, chi_value_AAA1, 'bo-', 'MarkerFaceColor', 'b','LineWidth',3.5);
plot(days_of_oxidation_AAA1,oxygen_content_AAA1,'bo-','MarkerFaceColor','b','Linewidth','3.5');

%xlabel('Days of oxidation');
xlabel('Oxygen Content (%)');
ylabel('\chi value');
%title('Chi Value vs Oxygen Content for AAA1 Asphalt System');
%text(oxygen_content_AAA1, chi_value_AAA1, num2str([chi_value_AAA1', oxygen_content_AAA1']), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
hold on

brown = [0.6 0.3 0];

set(gca,'LineWidth',2.5,'FontSize',16);
% Create the line plot for AAK1 system
%plot(days_of_oxidation_AAK1, chi_value_AAK1, 'go-', 'MarkerFaceColor', 'g','LineWidth',3.5);
%plot(oxygen_content_AAK1, chi_value_AAK1, 'go-', 'MarkerFaceColor', 'g','LineWidth',3.5);
plot(days_of_oxidation_AAK1,oxygen_content_AAK1,'go-','MarkerFaceColor','g','LineWidth',3.5);
%title('Chi Value vs Oxygen Content for AAA1 and AAK1 Asphalt Systems');

% Create the line plot for AAM1 system
%plot(days_of_oxidation_AAM1, chi_value_AAM1,'ro-', 'MarkerFaceColor', 'r','LineWidth',3.5);
%plot(days_of_oxidation_AAM1,oxygen_content_AAM1,'ro-','MarkerFaceColor','r');
plot(days_of_oxidation_AAM1,oxygen_content_AAM1,'ro-','MarkerFaceColor','r','LineWidth',3.5);
%plot(oxygen_content_AAM1, chi_value_AAM1, 'ro-', 'MarkerFaceColor', 'r','LineWidth',3.5);
%title('Chi Value vs Oxygen Content for All Three Asphalt Systems (Oxidized Asphaltenes)');

% Add legend
legend('AAA1', 'AAK1', 'AAM1');

% Create the box with table
%annotation('textbox', [0.65, 0.1, 0.25, 0.8], 'String', {'Asphalt System', 'Oxygen Content(%)', '/Chi Value', 'AAA1', oxygen_content_AAA1, chi_value_AAA1, 'AAK1', oxygen_content_AAK1, chi_value_AAK1, 'AAM1', oxygen_content_AAM1, chi_value_AAM1}, 'FitBoxToText', 'on');

% Adjust plot settings
xlim([0 30]);
ylim([0.10 0.3]);
grid off;
