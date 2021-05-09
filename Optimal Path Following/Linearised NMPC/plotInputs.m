% Plot history of inputs during simulation
figure
% Plot omega history
subplot(2,1,1)
stairs(times(1:length(times)), uHis(1,:))
title('Inputs')
ylabel("\omega [rad/s]")
grid on
% Plot a history
subplot(2,1,2)
stairs(times(1:length(times)), uHis(2,:))
grid on
ylabel("a [m/s^2]")
xlabel("Time (s)")
