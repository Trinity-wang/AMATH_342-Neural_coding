set(0,'defaultaxesfontsize',20);
set(0,'defaulttextfontsize',20); 


M=10000;

sample_list=rand(1,M);

[nlist,centerlist]=hist(sample_list,50); % nlist = number of samples per bin    centerlist = center of each bar in histogram

%just make a histogram of probabilities of counts lying in different bins

figure
bar(centerlist,nlist/(M)); % plots the center on x axis, + probabilty that a sample lands in each specific bin.


%Now interpret as a continuous probability distribution.  Why do we do this
%additional normalization (division)?

% deltac=centerlist(2)-centerlist(1);
% 
% figure
% bar(centerlist,nlist/(M*deltac));


