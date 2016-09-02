function [] = VEPcluster(AnimalName,Date)
% VEPcluster.m
%
%  Will take data from a retinotopic mapping experiment and try to cluster
%   / classify VEPs based on their waveforms
%INPUT: AnimalName - unique identifier for the animal as a number, e.g.
%            12345
%       Date - date of the experiment, e.g. 20160525

%OUTPUT: figures with scatter plots of VEP waveforms (maybe try PCA) based
%         on retinotopic position
%
% Created: 2016/08/31, Orange Line Metro, DC
%  Byron Price
% Updated: 2016/08/31
%  By: Byron Price

cd('~/CloudStation/ByronExp/Retino');

EphysFileName = sprintf('RetinoData%d_%d',Date,AnimalName); % no file identifier
    % because MyReadall does that for us
  
global centerVals Radius reps stimTime holdTime numStimuli w_pixels h_pixels ...
    DistToScreen numChans sampleFreq stimLen maxWin; %#ok<*REDEF>

StimulusFileName = sprintf('RetinoStim%d_%d.mat',Date,AnimalName);
load(StimulusFileName)

centerVals = stimParams.centerVals;
Radius = stimParams.Radius;
reps = stimParams.reps;
stimTime = stimParams.stimTime;
holdTime = stimParams.holdTime;
numStimuli = stimParams.numStimuli;
w_pixels = stimParams.w_pixels;
h_pixels = stimParams.h_pixels;
DistToScreen = stimParams.DistToScreen;

reps = reps-1;
% convert from allad to ChanData by filtering
[ChanData,timeStamps,tsevs,svStrobed] = ExtractSignal(EphysFileName);

% get LFP response to each stimulus (the VEPs)
[Response,meanResponse,~] = CollectVEPS(ChanData,timeStamps,tsevs,svStrobed);

[minVals,maxVals,minLatency,maxLatency] = GetStats(Response);

xPos = unique(centerVals(:,1));
yPos = unique(centerVals(:,2));
numX = length(xPos);
numY = length(yPos);
position = zeros(numStimuli,1);

mapping = zeros(numX,numY);
count = 1;
for ii=numY:-1:1
    for jj=1:numX
        mapping(jj,ii) = count;
        count = count+1;
    end
end
for ii=1:numStimuli
    xInd = find(logical(centerVals(ii,1) == xPos));
    yInd = find(logical(centerVals(ii,2) == yPos));
    position(ii) = mapping(xInd,yInd);
end

for ii=1:numChans
    w(1+2*(ii-1)) = figure();
    w(2+2*(ii-1)) = figure();
    for jj=1:numStimuli
%         Data = squeeze(Response(ii,jj,:,:));
%         [coeff,score] = pca(Data);
        figure(w(1+2*(ii-1)));
        Data = squeeze(minLatency(ii,jj,:));
        bootData = zeros(2000,2);
        for kk=1:2000
            indeces = random('Discrete Uniform',reps,[reps,1]);
            bootData(kk,1) = mad(Data(indeces),1);
            bootData(kk,2) = quantile(Data(indeces),0.75)-quantile(Data(indeces),0.25);
        end
        subplot(numY,numX,position(jj));histogram(bootData(:,1))
        %scatter(score(:,1),score(:,2))
%         scatter(squeeze(minLatency(ii,jj,:)),squeeze(minVals(ii,jj,:)));
        figure(w(2+2*(ii-1)));
        subplot(numY,numX,position(jj));histogram(bootData(:,2))
%         scatter(squeeze(minVals(ii,jj,:)),squeeze(minLatency(ii,jj,:)));
    end
end


end

function [ChanData,timeStamps,tsevs,svStrobed] = ExtractSignal(EphysFileName)
    % Extract LFP signals from allad, filter, get timestamps
    global numChans sampleFreq;
    % read in the .plx file

    if exist(strcat(EphysFileName,'.mat'),'file') ~= 2
        MyReadall(EphysFileName);
    end

    EphysFileName = strcat(EphysFileName,'.mat');
    load(EphysFileName)

    
    sampleFreq = adfreq;

    Chans = find(~cellfun(@isempty,allad));
    numChans = length(Chans);

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
end

function [Response,meanResponse,strobeTimes] = CollectVEPS(ChanData,timeStamps,tsevs,svStrobed)
    global numChans numStimuli reps stimLen stimTime sampleFreq;
    strobeStart = 33;
    strobeTimes = tsevs{1,strobeStart};
    stimLen = round(0.2*sampleFreq); % about 250 milliseconds
    smoothKernel = 4;
    % COLLECT DATA IN THE PRESENCE OF VISUAL STIMULI
    Response = zeros(numChans,numStimuli,reps,stimLen);
    meanResponse = zeros(numChans,numStimuli,stimLen);
    for ii=1:numChans
        for jj=1:numStimuli
            stimStrobes = strobeTimes(svStrobed == jj);
            for kk=1:reps
                stimOnset = stimStrobes(kk);
                [~,index] = min(abs(timeStamps-stimOnset));
                temp = ChanData(index:index+stimLen-1,ii);
                Response(ii,jj,kk,:) = temp;
            end
            meanResponse(ii,jj,:) = smooth(mean(squeeze(Response(ii,jj,:,:)),1),smoothKernel);
        end
    end

end

function [minVals,maxVals,minLatency,maxLatency] = GetStats(Response)
       
    global numChans numStimuli stimLen sampleFreq maxWin;
    maxWin = round(0.1*sampleFreq):stimLen;
    [minVals,minLatency] = min(Response,[],4);
    [maxVals,maxLatency] = max(Response(:,:,:,maxWin),[],4);
    
%     for ii=1:numChans
%         for jj=1:numStimuli
%             
%         end
%     end
    minLatency = minLatency./sampleFreq;
    maxLatency = (maxLatency+maxWin(1))./sampleFreq;
    
end

