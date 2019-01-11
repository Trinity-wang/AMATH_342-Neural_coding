%% Class Notes 1/8/18

% Neurons and Stuff Review
%     Input = electrical
%     turned into a voltage difference across membrane: action potentials
%     output = electrical

% What Kind of Data can we analyze?
%   Record spikes in membrane potential; measure spiking rate compared
%   against stimulus given; may allow us to see what this neuron is encoding,
%   i.e. orientation, color, etc for visual stimuli

    % Averaged Spike Statistics; Frequency of spikes
    % Precise Spike Timing: information carried in exact time spikes occur
    %                       after stimulus, pattern they form
    %
    
% Marr's 3 Levels
%   Computational Level (highest): what does system do? learning + behavior
%   algorithmic Level (middle): How does it do this? What math does it use?
%   Implementation (lowest): Specific representations of algorithms in
%                            neural hardware

% Creating a random variable
    % x = rand() generates a single sample of a random variable
        % pseudorandom, based on some seed number and a consistent formula
        % for us, NEED to develop our own random starting point
        rand('state', sum(100*clock));
        % Generate binary random variable with probability p;
        p = 0.25; % probability of a 1
        X = round(rand() + p - 0.5);
        % to make a list of N random variables rand(1, N);
        %Generate single spiketrain rand('state',sum(100*clock));
        nsec=1 ;
        T=1;
        deltat=0.001;
        r=100; % neuron firing rate
        p=r*deltat; % prob of a neuron firing in a given time window
        numbins=round(T/deltat);  % number of time intervals
        spiketrain=round(rand(1,numbins) + (p-1/2));
        figure; 
        imagesc(spiketrain) ...
            
    %Generate many "trials" of spiketrains
    numtrials=10; 
    spiketrain=round(rand(numtrials,numbins) + (p-1/2));
    figure; 
    imagesc(spiketrain) 
    xlabel('time') 
    ylabel('trial')
    
    %Compute the average spike rate, and standard deviation
    rate_per_trial=1/T * sum(spiketrain,2)
    mean_rate_per_trial = mean(rate_per_trial) 
    std_dev_rate_per_trial = std(rate_per_trial)
    