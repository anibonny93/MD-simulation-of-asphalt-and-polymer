% Load the data into a variable called 'data'
data = [chiall_1];

% Fit an ARIMA model to the data
model = arima(0,1,0); % (p,d,q) where p is the AR order, d is the degree of differencing, and q is the MA order
model = estimate(model, data);

% Forecast the next 500 data points
[y, yMSE] = forecast(model, 500);

% Plot the original data and the forecast
plot(1:length(data), data, 'o', 1:length(data)+500, [data y], '-');
legend('Original Data', 'Forecast');
xlabel('Time');
ylabel('Value');
