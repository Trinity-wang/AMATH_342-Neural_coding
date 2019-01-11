%% Matlab Tutorial

x = [1 2 3 4 5];
elements = x([2 1 3]);

%% 3.1

t = 1:10;
signal_vector = sin(t);
I_vector = cumsum(signal_vector);
plot(t, I_vector, '-*');
xlabel('time(s)');
ylabel('Signal Cummulative Sum');
title('time vs cummulative sum');

for i = 1:9
    disp(i);
end
%% 3.2
tspan = 0:30;
numberOn = zeros(1,length(tspan));
numberOn(1) = 1;
for t = 2:length(tspan);
    numberOn(t) = numberOn(t-1)*3;
end
plot(tspan, numberOn, 'r-');

%%  
diary on;
x = 1;
b = 1 + x;
x = 3 + 1

a = [-1 1 2 -2 3 -3 4 -4 5 -5];
aNonNeg = a(find(a >= 0));
a >= 0

%% Plotting In 3D space

syms s t
x = sin(s)*cos(t);
y = sin(s)*sin(t);
z = s;
ezsurf(x,y,z, [0, 4*pi, 0, 2*pi]);
%%
