% Plot history of positions during simulation as well as track
figure
plot(xiHis(1,:),xiHis(2,:))
hold on
% axis([track_xmin track_xmax track_ymin track_ymax])
plot(xCenter,yCenter,'r')
plot(xOuter,yOuter,'k')
plot(xInner,yInner,'k')
title('Positions')
xlabel('X')
ylabel('Y')
legend('vehicle','reference', 'boundary')