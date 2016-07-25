function [] = MapRetinotopyCallaway(AnimalName,Date,yesNo)
% MapRetinotopyCallaway.m
%
%  Will take data from a retinotopic mapping experiment and extract the
%   retinotopy of the LFP recording electrode. The stimulus used here is
%   the periodic drifting bar with flashing counter-phase checkerboard from
%   the file RetinotopyCallaway.m
%     Must have the Matlab Offline Files SDK on the current path
%
%INPUT: AnimalName - unique identifier for the animal as a number, e.g.
%            12345
%       Date - date of the experiment input as a number yearMonthDay, 
%            e.g. 20160525
%       Chans - channel numbers, input as [6,8], defaults to 6 and 8
%OUTPUT: a plot
%
% Created: 2016/05/31, 24 Cummington, Boston
%  Byron Price
% Updated: 2016/07/25
%  By: Byron Price

% read in the .plx file
EphysFileName = strcat('RetinoCallData',num2str(Date),'_',num2str(AnimalName));

if exist(strcat(EphysFileName,'.mat'),'file') ~= 2
    MyReadall(EphysFileName);
end

StimulusFileName = strcat('RetinoCallStim',num2str(Date),'_',num2str(AnimalName),'.mat');
EphysFileName = strcat(EphysFileName,'.mat');
load(EphysFileName)
load(StimulusFileName)

if nargin < 3
    yesNo = 1;
end

sampleFreq = adfreq;
% newFs = 5;
% downSampleRate = sampleFreq/newFs;


% tsevs are the strobed times of stimulus onset, then offset
%  Onset at tsevs{1,33}(2), offset at tsevs{1,33}(3), onset at
%  tsevs{1,33}(4), offset at 5, etc.
% allad contains the continuous data from each channel, which appear to be
%  recorded at 1000 Hz rather than 40,000

%totalAD = size(allad,2);
%totalSEVS = size(tsevs,2);

Chans = find(~cellfun(@isempty,allad));numChans = length(Chans);
strobeStart = 33;

% lowpass filter the data
dataLength = length(allad{1,Chans(1)});

ChanData = zeros(dataLength,numChans);
preAmpGain = 1/1000;
for ii=1:numChans
    voltage = ((allad{1,Chans(ii)}).*SlowPeakV)./(0.5*(2^SlowADResBits)*adgains(Chans(ii))*preAmpGain);
    temp = smooth(voltage,0.05*sampleFreq);
    n = 30;
    lowpass = 100/(sampleFreq/2); % fraction of Nyquist frequency
    blo = fir1(n,lowpass,'low',hamming(n+1));
    ChanData(:,ii) = filter(blo,1,temp);
end


timeStamps = 0:1/sampleFreq:dataLength/sampleFreq-1/sampleFreq;

if length(timeStamps) ~= dataLength
    display('Error: Review allad cell array and timing')
    return;
end
strobeData = tsevs{1,strobeStart};
stimStart = round(0.05*sampleFreq);
stimLen = round(0.05*sampleFreq);

% COLLECT DATA IN THE PRESENCE OF VISUAL STIMULI
numDirs = 4;
Response = cell(numChans,numDirs);
numFlashes = zeros(numDirs,1);

diffs = abs(diffs);
for ii=1:numChans
    for jj=1:numDirs
        numFlashes(jj) = sum(diffs);
        Response{ii,jj} = zeros(reps,numFlashes(jj),stimLen);
    end
end

for ii=1:numChans
    count = 1;
    for jj=1:numDirs
        for kk=1:reps
            for ll=1:numFlashes(jj)
                stimOnset = strobeData(count);
                [~,index] = min(abs(timeStamps-stimOnset));
                temp = ChanData(index+stimStart:index+stimStart+stimLen-1,ii);
                Response{ii,jj}(kk,ll,:) = temp;
                count = count+1;
            end
        end
    end
end

% BOOTSTRAP FOR 95% CONFIDENCE INTERVALS OF STATISTIC IN ABSENCE OF VISUAL STIMULI
%  AND STANDARD ERRORS
%  started trial with 120 seconds of a blank screen
N = 5000; % number of bootstrap samples
noStimLen = startPause*sampleFreq;

bootPrctile = zeros(numChans,1); % 99 percentile
bootStat = zeros(numChans,2);
for ii=1:numChans
    Tboot = zeros(N,1);
    for jj=1:N
        indeces = random('Discrete Uniform',noStimLen-(stimLen+stimStart)*2,[reps,1]);
        temp = zeros(reps,stimLen);
        for kk=1:reps
            temp(kk,:) = ChanData(indeces(kk)+stimStart:indeces(kk)+stimStart+stimLen-1,ii);
        end
        Tboot(jj) = abs(min(mean(temp,1)));
    end
    %figure();histogram(Tboot);
    bootPrctile(ii) = quantile(Tboot,1-1/100);
    bootStat(ii,1) = mean(Tboot);
    bootStat(ii,2) = std(Tboot);
end

dataStat = cell(numChans,numDirs);
for ii=1:numChans
    for jj=1:numDirs
        temp = squeeze(mean(Response{ii,jj},1));
        dataStat{ii,jj} = abs(min(temp,[],2));
    end
end

% BOOTSTRAP FOR STANDARD ERROR OF STATISTIC IN PRESENCE OF VISUAL STIMULI
N = 1000;
dataError = cell(numChans,numDirs);
for ii=1:numChans
    for jj=1:numDirs
        dataError{ii,jj} = zeros(numFlashes(jj),1);
        for kk=1:numFlashes(jj)
            Tboot = zeros(N,1);
            for ll=1:N
                indeces = random('Discrete Uniform',reps,[reps,1]);
                group = squeeze(Response{ii,jj}(indeces,kk,:));
                meanGroup = mean(group,1);
                Tboot(ll) = abs(min(meanGroup));
            end
            dataError{ii,jj}(kk) = std(Tboot);
        end
    end
end

% WALD TEST - VEP magnitude significantly greater in presence of a stimulus
%  than in the absence of a stimulus
significantStimuli = cell(numChans,numDirs);
alpha = 0.05/sum(numFlashes);
for ii=1:numChans
    for jj=1:numDirs
        significantStimuli{ii,jj} = zeros(numFlashes(jj),1);
        for kk=1:numFlashes(jj)
            W = (dataStat{ii,jj}(kk)-bootStat(ii,1))/sqrt(dataError{ii,jj}(kk)^2+bootStat(ii,2)^2);
            c = norminv(1-alpha,0,1);
            if W > c
                significantStimuli{ii,jj}(kk) = dataStat{ii,jj}(kk); % or equals W itself
            end
        end
    end    
end

if yesNo == 1
    figure();
end
count = 1;
for ii=1:numChans
    for jj=1:2
        x = centerPos{jj}(logical(diffs));y = centerPos{jj+2}(logical(diffs));
        stimVals = log(significantStimuli{ii,jj}*significantStimuli{ii,jj+2}');
        if yesNo == 1 && jj==1
            subplot(numChans,2,count);imagesc(x,y,stimVals);set(gca,'YDir','normal');h=colorbar;
            title(sprintf('Right/Up Sweep Retinotopy for Channel %d , Animal %d',ii,AnimalName));ylabel(h,'VEP Magnitude (\muV)');
            xlabel('Horizontal Screen Position (pixels)');ylabel('Vertical Screen Position (pixels)');
            colormap('jet');
        end
        if yesNo == 1 && jj==2
            subplot(numChans,2,count);imagesc(x,y,stimVals);set(gca,'YDir','normal');h=colorbar;
            title(sprintf('Left/Down Sweep Retinotopy for Channel %d , Animal %d',ii,AnimalName));ylabel(h,'VEP Magnitude (\muV)');
            xlabel('Horizontal Screen Position (pixels)');ylabel('Vertical Screen Position (pixels)');
            colormap('jet');
        end
        count = count+1;
    end
end
% if yesNo == 1
%     savefig(strcat('RetinoMap',num2str(Date),'_',num2str(AnimalName),'.fig'));
% end
% figure();
% hold on;
% xlim([0 w_pixels])
% ylim([0 h_pixels])
% for ii=1:numChans
%     plot(CenterPoints(1,ii),CenterPoints(2,ii),'*')
% end
end
