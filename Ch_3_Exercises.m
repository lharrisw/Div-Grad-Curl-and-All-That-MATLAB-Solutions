%% Example 1
clear
close all
clc

s = linspace(0,sqrt(2));
S = linspace(0,1);
x = @(s) s/sqrt(2);
y = @(s) s/sqrt(2);
X = @(s) s;
Y = @(s) s;

figure
hold on
grid on
plot(x(s),y(s),'linewidth',2)
plot(X(S),Y(S),'ro')

%% Example 2
clear
clc
close all

R = 3;
s = linspace(0,R);
x = @(s) s;
y = @(s) sqrt(R^2-s.^2);

figure
hold on
grid on
plot(x(s),y(s))
axis equal

%% Example 4
clear
clc
close all

s = linspace(-1,1);
x = s;
yplus = sqrt(1-s.^2);
yminus = -sqrt(1-s.^2);

figure
hold on
grid on
plot(x,yplus)
plot(x,yminus)