clc;
close all;
clear all;

phi = [0.001:0.001:0.999, 0.9999, 0.99999, 0.999999, 0.99999999, 0.999999999, 0.9999999999999999];
k = 0.001985875;
chi = 0.1223;
N1 = 10;
N2 = 10;
X = [];
Y = [];

T_values = 100:10:800;

for T = T_values
    F = @(phi) T.*(((phi .* log(phi)) + ((1 - phi) ./ N2) .* log(1 - phi) + chi .* phi .* (1 - phi)));
    y = F(phi);

    k = convhull(phi, y);
    [phik, yk] = deal(phi(k), y(k));
    j = diff(phik) > 0;
    [phik, yk] = deal(phik(j), yk(j));
    interDist = vecnorm(diff([phik(:), yk(:)], 1, 1), 2, 2);
    [~, imax] = max(interDist);
    x1 = phik(imax);
    y1 = yk(imax); %the two tangent points (x1,y1) and (x2,y2)
    x2 = phik(imax + 1);
    y2 = yk(imax + 1);

    X = [X, x1, x2];
    Y = [Y, T, T];
end

% Truncate X to match the length of Y
%X = X(1:length(Y));

plot(X, Y, 'o');
xlabel('X');
ylabel('Y');
title('Plot of X versus Y');

