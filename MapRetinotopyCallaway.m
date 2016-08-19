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
% Updated: 2016/08/18
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

Chans = find(~cellfun(@isempty,allad));numChans = length(Chans);
strobeStart = 33;

% lowpass filter the data
dataLength = length(allad{1,Chans(1)});

ChanData = zeros(dataLength,numChans);
preAmpGain = 1;
for ii=1:numChans
    voltage = 1000.*((allad{1,Chans(ii)}).*SlowPeakV)./(0.5*(2^SlowADResBits)*adgains(Chans(ii))*preAmpGain);
    n = 30;
    lowpass = 100/(sampleFreq/2); % fraction of Nyquist frequency
    blo = fir1(n,lowpass,'low',hamming(n+1));
    ChanData(:,ii) = filter(blo,1,voltage);
end


timeStamps = 0:1/sampleFreq:dataLength/sampleFreq-1/sampleFreq;

if length(timeStamps) ~= dataLength
    display('Error: Review allad cell array and timing')
    return;
end
strobeTimes = tsevs{1,strobeStart};
stimStart = round(0*sampleFreq);
stimLen = round(0.2*sampleFreq);
minWin = round(0.04*sampleFreq):1:round(.1*sampleFreq);
maxWin = round(0.1*sampleFreq):1:round(.2*sampleFreq);

% COLLECT DATA IN THE PRESENCE OF VISUAL STIMULI
Response = cell(numChans,numDirs);
meanResponse = cell(numChans,numDirs);
numFlashes = zeros(numDirs,1);
smoothKernel = 4;

for ii=1:numChans
    for jj=1:numDirs
        numFlashes(jj) = sum(svStrobed == jj);
        Response{ii,jj} = zeros(numFlashes(jj),reps,stimLen);
        meanResponse{ii,jj} = zeros(numFlashes(jj),stimLen);
    end
end

for ii=1:numChans
    for jj=1:numDirs
        dirStrobes = strobeTimes(svStrobed==jj);
        count = 1;
        for kk=1:numFlashes(jj)
            for ll=1:reps
                stimOnset = dirStrobes(count);
                [~,index] = min(abs(timeStamps-stimOnset));
                temp = ChanData(index+stimStart:index+stimStart+stimLen-1,ii);
                Response{ii,jj}(kk,ll,:) = temp;
                count = count+1;
                clear temp;
            end
            temp = mean(squeeze(Response{ii,jj}(kk,:,:)),1);
            meanResponse{ii,jj}(kk,:) = smooth(temp,smoothKernel);
        end
    end
end

% BOOTSTRAP FOR 95% CONFIDENCE INTERVALS OF STATISTIC IN ABSENCE OF VISUAL STIMULI
%  AND STANDARD ERRORS
%  intermixed trial with 30 seconds of grey screen
N = 1000; % number of bootstrap samples
noStimLen = holdTime*sampleFreq-(stimLen+stimStart)*2;

baseStats(1:numChans) = struct;
for ii=1:numChans
    Tboot = zeros(N,1);
    for jj=1:N
        pauseOnset = strobeTimes(svStrobed == 0);
        val = random('Discrete Uniform',length(pauseOnset),1);
        [~,index] = min(abs(timeStamps-pauseOnset(val)));
        indeces = random('Discrete Uniform',noStimLen,[reps,1]);
        temp = zeros(reps,stimLen);
        for kk=1:reps
            temp(kk,:) = ChanData(index+indeces(kk)+stimStart:index+indeces(kk)+stimStart+stimLen-1,ii);
        end
        Tboot(jj) = abs(min(mean(temp(:,minWin),1)));
    end
    %figure();histogram(Tboot);
    baseStats(ii).lbound = quantile(Tboot,1/100);
    baseStats(ii).ubound = quantile(Tboot,1-1/100);
    baseStats(ii).mean = mean(Tboot);
    baseStats(ii).sem = std(Tboot);
end

dataStats(1:numChans) = struct;
for ii=1:numChans
    dataStats(ii).mean = cell(numDirs,1);
    dataStats(ii).sem = cell(numDirs,1);
    for jj=1:numDirs
        dataStats(ii).mean{jj} = zeros(numFlashes(jj),1);
        dataStats(ii).sem{jj} = zeros(numFlashes(jj),1);
        for kk=1:numFlashes(jj)
            temp = squeeze(meanResponse{ii,jj}(kk,:));
            dataStats(ii).mean{jj}(kk) = abs(min(temp(minWin)));
        end
    end
end

% BOOTSTRAP FOR STANDARD ERROR OF STATISTIC IN PRESENCE OF VISUAL STIMULI
N = 1000;
for ii=1:numChans
    for jj=1:numDirs
        dataStats(ii).sem{jj} = zeros(numFlashes(jj),1);
        for kk=1:numFlashes(jj)
            Tboot = zeros(N,1);
            for ll=1:N
                indeces = random('Discrete Uniform',reps,[reps,1]);
                group = squeeze(Response{ii,jj}(indeces,kk,:));
                meanGroup = mean(group,1);
                Tboot(ll) = abs(min(meanGroup(minWin)));
            end
            dataStats(ii).sem{jj}(kk) = std(Tboot);
        end
    end
end

% WALD TEST - VEP magnitude significantly greater in presence of a stimulus
%  than in the absence of a stimulus
significantStimuli = cell(numChans,numDirs);
alpha = 0.05/sum(numFlashes);
c = norminv(1-alpha,0,1);
for ii=1:numChans
    for jj=1:numDirs
        significantStimuli{ii,jj} = zeros(numFlashes(jj),1);
        for kk=1:numFlashes(jj)
            W = (dataStats(ii).mean{jj}(kk)-baseStats(ii).mean)/...
                sqrt(dataStats(ii).sem{jj}(kk)^2+baseStats(ii).sem^2);
            if W > c
                significantStimuli{ii,jj}(kk) = dataStats(ii).mean{jj}(kk); % or equals W itself
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
