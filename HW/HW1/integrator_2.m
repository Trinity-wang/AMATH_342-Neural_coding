tspan = 0:0.01:10;
signal_vector = sin(tspan);
thresh = 100;
sum = 0;
for index = 1:length(tspan)
    if sum > thresh
        break
    else
        sum = sum + signal_vector(index);
    end
end
timeOverTresh = tspan(index)