clear all;

load x
load y

axis equal

nChar = length(x) - 1;
plot(x , y);
title('Number of Characteristics: ', nChar);
xlabel('X/Throat');
ylabel('Y/Throat');
xlim([0 5]);
ylim([0 2]);



%z = (y(51) - y(50))/(x(51) - x(50));
