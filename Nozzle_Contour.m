clear all;

load x
load y


nChar = length(x) - 1;
axis equal
plot(x , y);
axis equal

title('Number of Characteristics: ', nChar);
xlabel('X/Throat');
ylabel('Y/Throat');




%z = (y(51) - y(50))/(x(51) - x(50));
