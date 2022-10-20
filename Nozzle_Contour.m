clear all;

load x
load y

axis equal

plot(x , y);
xlim([0 40]);
ylim([0 40]);

z = (y(51) - y(50))/(x(51) - x(50));
