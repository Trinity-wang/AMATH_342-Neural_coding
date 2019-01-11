% Exercise 3.2 
tspan = 0:30;
number_on = zeros(length(tspan),1);
number_on(1)  = 1;
for t = 1:length(tspan) - 1
    number_on(t + 1) = 3 * number_on(t);
end
plot(tspan, number_on, 'ro-');
xlabel('time(s)');
ylabel('neurons on');

