%% HW2 Visual Stimuli


%% Part 1: Finding a good STA

close all; clear all; clc;

exp_sec = 1000; % length of experiment in seconds
run('v1_white_noise_exp.m')

% stim will be a 11 X 11 X second*FrameRate frames array where each of the Frames is one image shown to our V1 neuron

% spikeTrain will be a vector 1 X seconds*frameRate long where a 1
% represents a spike at that specific time and a 0 = no spike.

% Part 1.1
frameRate = 60; % frames per second.
spikeIndexes = (find(spikeTrain))';

% delayTime should be between 0 and 0.5 seconds
for j = 0:20
delayTime = 0.025*j; % seconds to look back
lookback = ceil(frameRate*delayTime); % frames to lookBack

stimsum = zeros(11,11);
STA = zeros(11,11);

for index = spikeIndexes
    stimsum = stimsum + stim(:,:, index - lookback);
end
STA = stimsum/length(spikeIndexes);

subplot(3,7, j + 1);
colormap('gray');
imagesc(STA);
title(compose('%.3f seconds delay', delayTime), 'Fontsize', 12);
end

%% Part 2
clear all; close all; clc;

% Part 2.1

% This part of the code generates the firing rate distributions for the
% stimulus values of interest.

ntrials = 10000;
cell_num = 1;
midStim = 20; % What stimulus our analysis is centered around
dx = -5; % How far below/above our center we use as a secondary stimulus value for comparison


stimDir = midStim;

run('generate_noisy_data_cockroach.m');

binwidth = 1; % binwidth for all the histograms


spikeCounts1 = sum(spiketrain'); % total count of spikes in each trial
histogram(spikeCounts1, 'BinWidth', binwidth, 'Normalization', 'probability', 'FaceColor', 'red');

hold on;


edges = 0:binwidth:150;
trialCounts = zeros(2, length(edges) - 1);

% stores the number of trials that occured with this given stimulus within
% each of the bins whose edges are specified in 'edges'
trialCounts(1,:) = histcounts(spikeCounts1, edges); 



stimDir = midStim + dx;

run('generate_noisy_data_cockroach.m');

spikeCounts2 = sum(spiketrain'); % total count of spikes in each trial
histogram(spikeCounts2, 'BinWidth', binwidth, 'Normalization', 'probability', 'FaceColor', 'blue');

hold on;

trialCounts(2,:) = histcounts(spikeCounts2, edges);


set(gca,'fontsize',20)
name1 = compose('%d degrees', midStim);
name2 = compose('%d degrees', midStim + dx);
legend([name1 name2], 'Fontsize', 25);
title('Total Spike Count Distributions For Varying Stimulus Values', 'Fontsize', 25);
xlabel('Total Spikes in a Trial', 'Fontsize', 20);
ylabel('Proportion of Trials', 'Fontsize', 20);
hold off;



% This part of the code analyzes the distributions and, using maximum likelihood discrimination, attempts to
% guess the stimulus that created a given firing rate, also providing the
% error rate of such a guess.

% Part 2.1a)

% finds good threshold value
midIsGreater = (trialCounts(1,:) >= trialCounts(2,:)).*(trialCounts(1,:) > 5);
otherIsGreater = (trialCounts(2,:) > trialCounts(1,:)).*(trialCounts(2,:) > 5);
[index dxIsLess] = max([find(midIsGreater > 0, 1), find(otherIsGreater > 0,1)]);
%dxIsLess = 1 if mid + dx dist comes first and == 2 if mid comes first
threshold = edges(index); % threshold firing rate used in determining which guess to make.
   


    
    % Find error rate
    errorRate = 0;
    Pmid = 0.5; % chance of the stimulus being mid degrees
    PmidPlusX = 0.5; % chance of stimulus being mid + dx degrees
    
    if dxIsLess == 1
       % error due to guessing 50 + dx when right guess was 50
       errorRate = errorRate + sum(trialCounts(1, 1:index - 1)) / ntrials * (Pmid); 
       % error due to guessing 50 when right guess was 50 + x
       errorRate = errorRate + sum(trialCounts(2, index:end)) / ntrials * (PmidPlusX);   
    else
       % error due to guessing 50 + dx when right guess was 50
       errorRate = errorRate + sum(trialCounts(1, index:end)) / ntrials * (Pmid); 
       % error due to guessing 50 when right guess was 50 + x
       errorRate = errorRate + sum(trialCounts(2, 1:index - 1)) / ntrials * (PmidPlusX); 
    end
    
    threshold
    errorRate
    fractionCorrect = 1 - errorRate
    
    midIsLess = find(trialCounts(1,:) <= trialCounts(2,:));
    midIsMore = find(trialCounts(2, :) < trialCounts(1,:));
    errorRate2 = sum(trialCounts(1,midIsLess)/ntrials * Pmid) + sum(trialCounts(2,midIsMore) / ntrials * PmidPlusX)
    
% Part 2.1b)
    % With dx = -5 we can get a 15% error rate, but any lower than this and
    % the error starts going up again because the tuning curve is symmetric
    % around 45 degrees as we saw in hw 1, so dx must be MUCH lower to get error < 10%
    % somewhere between -13 and -14 looks right
    
    % With dx between 3 and 4 we get ~10 % error
    
    
% Part 2.1c)
    % Now set midStim = 20
    
    % With dx < 0, there will never be a point where there will be a
    % noticeable difference because, as we saw in HW 1, the tuning curve
    % for Cell 1 flat lines at anything 30 degrees and below.
    
    % With dx between 15 and 16, it works well. Error is ~ 10%

    
% Part 2.2

    % change mean firing rate so it is 20...