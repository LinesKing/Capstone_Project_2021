% Plot history of states during simulation
% Plot x history
figure
subplot(4,1,1)
plot(times(1:length(times)), xiHis(1,1:end-1))
title('States')
ylabel("x [m]")
grid on
% Plot y history
subplot(4,1,2)
plot(times(1:length(times)), xiHis(2,1:end-1))
ylabel("y [m]")
grid on
% Plot v history
subplot(4,1,3)
plot(times(1:length(times)), xiHis(3,1:end-1))
ylabel("\theta (\circ)")
grid on
% Plot theta history
subplot(4,1,4)
plot(times(1:length(times)), xiHis(4,1:end-1))
grid on
ylabel("v [m/s]")
