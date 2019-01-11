% 
% This code will generate a spiking response to an input based on a tuning
% curve and assuming (inhomogeneous) Poisson firing. The response adapts
% over time, with a time constant of tau
%

clear all; close all; clc;

rand('state',sum(100*clock));
time1 = clock;

% stimDir = input('Input the direction of your stimulus in degrees ');
% cell_num = input('Which cell do you want to record from (1,2,3) ');
% ntrials = input('How many repeated trials would you like to perform? ');

% MY CODE
ntrials = 300;  % number of trials to run each stimulus at
nmsec = 300; % number of milliseconds to record for
stims = 0:10:90;
numStims = length(stims);   % number of bins the 0-90 degrees stimulus is separated into

for m = 1:3
    cell_num = m;
    trialAveRate = zeros(numStims, nmsec);
    meanRates = zeros(numStims,ntrials);
for k = 1:numStims
    stimDir = stims(k);


% NOT MY CODE
times= 1:nmsec; % vector of time points (1 msec apart)

spiketrain = zeros(ntrials,nmsec);      % set up output data

rate = cockroach_tuning(stimDir, cell_num); %returns rate, in Hz.       
tau = 100;      % adaptation time constant in msec
delta_t=0.001; %time bin, in seconds (1 msec)
ratefun = rate*exp(-times/tau);  % adapting rate function 

spiketrain = round(rand(ntrials,nmsec) + ones(ntrials,1)*((ratefun*delta_t)) - 1/2);

% MY CODE
figure(1); %Plots raster plots for each of given stimulus values across all 3 cell types
subplot(3,numStims, numStims*(m - 1) + k);
imagesc(spiketrain);
str = compose('Type %d, %d%s', [cell_num stimDir char(176)]);
title(str);
xlabel('time(ms)');
ylabel('trial number');

% average firing rate per second across trials at each individual timepoint
 trialAveRate(k,:) = sum(spiketrain)/(nmsec * delta_t);
 % average firing rate per second for each trial across entire time period
 meanRates(k, :) = mean(spiketrain')/delta_t;

figure(2); % Plots trial average firing rate over time to see how it fluctuates/trends
subplot(3, numStims, numStims*(m - 1) + k);
plot(times, trialAveRate(k,:));
title(str);
ylim([0 350]);
xlim([0 300])
xlabel('time(ms)');
ylabel('trial averaged firing rate (per sec)');
end

% Plots average (across entire time interval) firing rate vs stimulus angle with error bars
% Error measures standard deviation between mean firing rate (across entire
% interval) for each of the many trials.
figure(3); 
subplot(3,1, m);
x  = stims; % All the stimulus values (angles) the cell was tested at
means = mean(meanRates.'); % average firing rate over entire time interval
error = std(meanRates.');
colors = ['r','b','k'];
errorbar(x,means, error, colors(m));
title(compose('Cell Type %d', m));

xlabel(compose('Stimulus Angle (%s)', char(176)));
ylabel('Mean Firing Rate (per sec)');
ylim([0,200]);
xlim([0,90]);
end