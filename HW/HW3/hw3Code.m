% AMATH 342 HW #3

%%  Part 1: Filtering of Inputs: what matters in driving the membrane response?
clear all; close all; clc;

% 1.1

deltat=0.2 ; %timestep
Tmax=35;

tlist=linspace(0,Tmax,Tmax/deltat + 1) ;

R = 2;
C = 5;

Iapplist1 = ones(1,length(tlist))*5.25372;    % I1(t)
Iapplist2= zeros(1,length(tlist));    % I2(t)
    for j = 1:length(tlist)
        if mod(j - 1,5) == 0
           Iapplist2(j) = 5.2084*5.25364; 
        end
    end   
    
run('euler_illustrateRC_two_inputs');
% Creates two vectors, Vlist1 and Vlist2, that give the voltage values over
% time in response to the respective currents.

tNow = 30;
V1at30  = Vlist1(find(tlist == tNow))
V2at30 = Vlist2(find(tlist == tNow))

% 1.2

% The explict solution for V(t) at a given time represents the sum of past
% values of I(t), weighted by a specific weighting function K(t'). Thus
% individual I(t) values do not matter on their own, but only the summed
% and weighted values of I(t) matter in determining voltage.

%% Part 2: Summation of Simulataneous Impulses: Sum linearly, sublinearly ,or superlinearly?

% 2.1
% Current-drive (RC) circuit
clear all; close all; clc;

thresh = 10;

magn = 5; % Impulse magnitude (microAmps?);
impulseWidth = 5; % Impulse width; milliseconds (must be divisible by deltat)

deltat = 0.1; % timestep; milliseconds

Tmax=10;
V0 = 0;

tlist=linspace(0,Tmax,Tmax/deltat + 1) ;

% RC = 10 ms
R = 2;
C = 5;

Iapplist = zeros(1,length(tlist)); % I(t)
    start = 1/deltat + 1;
    Iapplist(start:  start + impulseWidth/deltat) = magn;
    
run('euler_illustrateRC'); % modified version for this code.
% Creates a vector Vlist that gives voltage values over time in response to
% the given current, Iapplist.

% determining N

peakVolt = max(Vlist)
f = peakVolt / thresh
N = 1;
fNew = f;
names = ['N = 1'];

figure(1)
while fNew < 1
   N = N + 1
   names = [names compose('N = %d', N)];
   Iapplist = zeros(1, length(tlist));
   newImpulse = N*magn;
   Iapplist(start : start + impulseWidth/deltat) = newImpulse;
    
    run('euler_illustrateRC');
   
    
    
    peakVolt = max(Vlist);
    fNew = peakVolt / thresh
end
subplot(211);
plot(tlist, ones(1,length(tlist))*thresh, 'g');
legend(names, 'threshold');

subplot(212);
legend(names);



fNew
N
f

% N is the lowest value to reach above threshold...
% N and f are related by the fact that N is the smallest integer such that
% f * N > 1. We see for each trial that fNew = f * N.

% Use the explicit formula for for V(t) for an RC circuit to solve for N...
% int (kt It)  = peakVolt      (1)
% int (kt N*It) >= threshold
% N* int(kt It) >= theshold
% Use equation (1)
% N >= threshold / peakVolt = 1/f



%%
clear all; close all; clc;

% 2.2
% Conductance-based Input Model

%circuit parameters
R=2;
C=5;
E=11;
V0=0;

deltat=0.1 ; %timestep
Tmax=10;

tlist=linspace(0,Tmax,Tmax/deltat +1); %not my code

% define input conductance
magn = 0.5;
duration = 2; 

thresh = 10; % threshold voltage value

gapplist=zeros(1,length(tlist)); % our new "I(t)"
    start = 1/deltat + 1;
    gapplist(start :  start + duration/deltat) = magn;
    
    
run('euler_illustrateRC_conduct');
% Creates a vector Vlist of the voltage values, V(t), at each given timepoint in
% tlist.

peakVolt = max(Vlist)
f = peakVolt / thresh
N = 1;
fNew = f;
names = ['N = 1'];

figure (1)
while fNew < 1
   N = N + 1
   names = [names compose('N = %d', N)];
      
   gapplist = zeros(1, length(tlist));
   newImpulse = N*magn;
   gapplist(start:  start + duration/deltat) = newImpulse;
   
    run('euler_illustrateRC_conduct');
    peakVolt = max(Vlist);
    fNew = peakVolt / thresh
end
subplot(211);
plot(tlist, ones(1,length(tlist))*thresh, 'g');
legend(names, 'threshold');

subplot(212);
legend(names);

% these g(t) values do not sum linearly. Show this with the various values
% of N and fNew
% Before, N was related to f by the simple formula that N was the smallest
% integer value that satisfies N*f > 1, but now that is not the case.


%Looking at v(t) = integral ... for conductance based circuit we see it is a nonlinear equation...

%% Part 3: Hodgkin-Huxley Model of a Neuron

clear all; close all; clc;

% 3.1: firing rate with a constant current

thresh = -55; % threshold membrane potential before an spike is triggered
npoints=50000;
dt=0.01;        %timestep between each point in seconds 

Imin = 0;
Imax = 50;
stepSize = 0.05;

fireRates = zeros((Imax - Imin)/stepSize + 1, 1);

index = 1;
for k = Imin/stepSize : Imax / stepSize
 
I = ones(npoints,1) * k * stepSize; % constant current of given magn

run('HH'); 
% creates a vector v which represents the voltage of the cell in
% response to the given current, I, over time.

% calculate firing rate (per second)
    count = 0; % will count number of firing events in the given time interval
    for j = 1:length(v) - 1
       if v(j) < thresh && v(j + 1) > thresh % counts points where v jumps above threshold
           count = count + 1;
       end
    end

    fireRates(index) = (count /(npoints*dt))*1000; 
    % npoints * dt represents total amount of time that passes during trial
    % I multipy by 1000 to convert ms to seconds...
    index = index + 1;
end

% plot tuning curve I vs firing rate
figure
plot([Imin:stepSize:Imax], fireRates, 'r.-');
xlabel('Current I(t)');
ylabel('Firing rate (per second)');
ylim([0 200]);

%% 3.2: firing rate with sinusoidal current

clear all; close all; clc;

thresh = -55; % threshold membrane potential before an spike is triggered
npoints=50000;
dt=0.01;        %timestep between each point in seconds
tlist = [0:dt:(npoints - 1)*dt]'; 

Imin = 5.5;
Imax = 7;
stepSize = 0.001;

fireRates = zeros((Imax - Imin)/stepSize + 1, 1);
background = zeros(npoints, 1);

% define constants for sinusoidal background current
epsilon = 4;
omega = 2*pi/0.5;
background = epsilon*sin(2*pi*omega*tlist); % sinusoidal background current

index = 1;
for k = Imin/stepSize : Imax / stepSize
 
Imagn = k*stepSize;      
I = ones(npoints,1) * Imagn; % constant current of given magn
I = I + background; % sinusoidal background current

run('HH'); % creates a vector v which represents the voltage of the cell in response to the given current, I, over time.

% calculate firing rate (per second)
    count = 0; % will count number of firing events in the given time interval
    for j = 1:length(v) - 1
       if v(j) < thresh && v(j + 1) > thresh % counts points where v jumps above threshold
           count = count + 1;
       end
    end

    fireRates(index) = count / (npoints*dt)*1000; 
    % npoints * dt represents total amount of time that passes during trial
    % I multipy by 1000 to convert ms to seconds...
    index = index + 1;
end

%% plot tuning curve I = 0 : max
figure
plot([Imin:stepSize:Imax], fireRates, 'r');
ylim([0 max(fireRates + 15)]);
xlabel('Applied Current Magnitude (uA)');
ylabel('Firing Rate (per second)');

% near the threshold point, the background current has an interesting
% effect where it can either resonate with the increasing voltage, giving it
% just enough current to reach threshold, or it can counteract a rising
% voltage with a lower current and cut short the spike.

%% 3.3: add noise?

clear all; close all; clc;


trials = 20;

thresh = -55; % threshold membrane potential before an spike is triggered
npoints=10000;
dt=0.05;        %timestep between each point in seconds
tlist = [0:dt:(npoints - 1)*dt]'; 

Imin = 0;
Imax = 20;
stepSize = 0.1;
noiseMaxMagn = 0.05;

spikeCounts = zeros((Imax - Imin)/stepSize + 1, trials);

for trialNum = 1:trials

index = 1;
for k = Imin/stepSize : Imax / stepSize

Imagn = k*stepSize; 


I = zeros(npoints,1);
I(1) = Imagn;
rand('state', sum(clock*100));

for i = 1:length(I) - 1
   I(i + 1) = I(i) + 2*(rand() - 0.5)*noiseMaxMagn; 
end

run('HH');
% creates a vector v which represents the voltage of the cell in
% response to the given current, I, over time.

% calculate firing rate (per second)
    count = 0; % will count number of firing events in the given time interval
    for j = 1:length(v) - 1
       if v(j) < thresh && v(j + 1) > thresh % counts points where v jumps above threshold
           count = count + 1;
       end
    end

    spikeCounts(index, trialNum) = count; 
    % npoints * dt represents total amount of time that passes during trial
    % I multipy by 1000 to convert ms to seconds...
    index = index + 1;
end
end

fireRates = spikeCounts/ (npoints*dt)*1000;
% Finding fano factor: This is only relevant in areas where the added noise
% is enough to bump the current over the threshold required to trigger a
% spike. If the applied current is too low, there will always be no spikes
% even with the noise leading to zero variance amongst trials (giving a zero fano factor) and with
% too much applied current, there will ALWAYS be spikes, leading again to
% no variance between trials (giving a zero fano factor).

means = mean(spikeCounts');
vars = var(spikeCounts');

bestFit = polyfit(means, vars, 1);
bestVals = polyval(bestFit, means);
fano = bestFit(1)
%%
figure
plot(means, vars, 'ro'), hold on;
% plot(means, bestVals, 'b');
title('Determining the Fano Factor');
xlabel('Mean spike count');
ylabel('Variance in spike count');

%% For a single I, plot fano factor vs noise magn

clear all; close all; clc;

trials = 100;

thresh = -55; % threshold membrane potential before an spike is triggered
npoints=10000;
dt=0.05;        %timestep between each point in seconds

noiseMagnMin = 0;
noiseMagnMax = 1;
stepSize = 0.01;

spikeCounts = zeros((noiseMagnMax - noiseMagnMin)/stepSize, trials);
fanos = zeros((noiseMagnMax - noiseMagnMin)/stepSize, 1);
Imagn = 4;

index = 1;
for noiseMagn = noiseMagnMin : stepSize : noiseMagnMax
for trialNum = 1:trials

I = zeros(npoints,1);
I(1) = Imagn;
rand('state', sum(clock*100));
for i = 1:length(I) - 1
   I(i + 1) = I(i) + 2*(rand() - 0.5)*noiseMagn; 
end

run('HH');
% calculate firing rate (per second)
    count = 0; % will count number of firing events in the given time interval
    for j = 1:length(v) - 1
       if v(j) < thresh && v(j + 1) > thresh % counts points where v jumps above threshold
           count = count + 1;
       end
    end
    
    spikeCounts(index, trialNum) = count; 
    % npoints * dt represents total amount of time that passes during trial
    % I multipy by 1000 to convert ms to seconds...
    
end
    fanos(index) = var(spikeCounts(index, :)) / mean(spikeCounts(index, :));
    index = index + 1;
end
%
plot([noiseMagnMin : stepSize : noiseMagnMax], fanos, 'r');
ylabel('Fano Factor');
xlabel('noise magnitude');
title(compose('I(t) magnitude = %d', Imagn));