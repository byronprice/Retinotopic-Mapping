function [] = MapRetinotopyCallaway(AnimalName,Date)
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
%OUTPUT: plots
%
% Created: 2016/05/31, 24 Cummington, Boston
%  Byron Price
% Updated: 2017/10/22
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

driftSpeed = stimParams.driftSpeed; % units of pixels per second
% stimFreq = stimParams.stimFreq;
% Width = stimParams.Width;
w_pixels = stimParams.w_pixels;
h_pixels = stimParams.h_pixels;
reps = stimParams.reps;
checkRefresh = stimParams.checkRefresh;
holdTime = stimParams.holdTime;
driftTime = stimParams.driftTime;
centerPos = stimParams.centerPos; 
Flashes = stimParams.Flashes;
numDirs = stimParams.numDirs;
DirNames = stimParams.DirNames;
% DistToScreen = stimParams.DistToScreen;

stimulationFrequency = 1/checkRefresh;

sampleFreq = adfreq;

Chans = find(~cellfun(@isempty,allad));numChans = length(Chans);
strobeStart = 33;

% lowpass filter the data
dataLength = length(allad{1,Chans(1)});

ChanData = zeros(dataLength,numChans);
preAmpGain = 1;
for ii=1:numChans
    voltage = 1000.*((allad{1,Chans(ii)}).*SlowPeakV)./(0.5*(2^SlowADResBits)*adgains(Chans(ii))*preAmpGain);
    n = 10;
    lowpass = 100/(sampleFreq/2); % fraction of Nyquist frequency
    blo = fir1(n,lowpass,'low',hamming(n+1));
    temp = filter(blo,1,voltage);
    
    notch = 60/(sampleFreq/2);
    bw = notch/n;
    [b,a] = iirnotch(notch,bw);
    ChanData(:,ii) = filter(b,a,temp);
end


timeStamps = 0:1/sampleFreq:dataLength/sampleFreq-1/sampleFreq;

if length(timeStamps) ~= dataLength
    display('Error: Review allad cell array and timing')
    return;
end
strobeTimes = tsevs{1,strobeStart};
stimLen = round(driftTime*sampleFreq);

% COLLECT DATA IN THE PRESENCE OF VISUAL STIMULI
Response = cell(numChans,3);

for ii=1:numChans
    Response{ii,1} = zeros(2*reps,stimLen);
    Response{ii,2} = zeros(2*reps,stimLen);
    
    greyStrobes = strobeTimes(svStrobed==0);
    numGrey = length(greyStrobes);
    Response{ii,3} = zeros(2*numGrey,stimLen);
    
    count = 1;
    for jj=1:numGrey
       onsetTime = round(greyStrobes(jj)*sampleFreq); 
       offsetTime = onsetTime+stimLen-1;
       Response{ii,3}(count,:) = ChanData(onsetTime:offsetTime,ii);
       count = count+1;
       onsetTime = offsetTime+1;
       offsetTime = onsetTime+stimLen-1;
       Response{ii,3}(count,:) = ChanData(onsetTime:offsetTime,ii);
       count = count+1;
    end
    
    horzCount = 1;
    vertCount = 1;
    for jj=1:numDirs
        strobeNum = jj;
        currentStrobeTimes = strobeTimes(svStrobed==strobeNum);
        for kk=1:reps
           onsetTime = round(currentStrobeTimes(kk)*sampleFreq);
           offsetTime = onsetTime+stimLen-1;
           if strcmp(DirNames{jj},'Left') == 1 || strcmp(DirNames{jj},'Down') == 1
              tempLFP = flipud(ChanData(onsetTime:offsetTime,ii));
           else
              tempLFP = ChanData(onsetTime:offsetTime,ii);
           end
           
           if strcmp(DirNames{jj},'Left') == 1 || strcmp(DirNames{jj},'Right') == 1
                Response{ii,1}(horzCount,:) = tempLFP;
                horzCount = horzCount+1;
           elseif strcmp(DirNames{jj},'Down') == 1 || strcmp(DirNames{jj},'Up') == 1
                Response{ii,2}(vertCount,:) = tempLFP;
                vertCount = vertCount+1;
           end
        end
    end
end

time = linspace(0,driftTime,stimLen);
horzPosition = linspace(0,w_pixels,stimLen); %w_pixels
vertPosition = linspace(0,h_pixels,stimLen); %h_pixels

waveletSize = 8;
x = linspace(-waveletSize/2*checkRefresh,waveletSize/2*checkRefresh,round(waveletSize*checkRefresh*sampleFreq));
kernel = exp(-2*pi*x*1i*stimulationFrequency);
transformBaseline = zeros(numChans,1);
for ii=1:numChans
    numGrey = size(Response{ii,3},1);
    temp = zeros(numGrey,stimLen);
    for jj=1:numGrey
        data = Response{ii,3}(jj,:);
        data = conv(data,kernel,'same');
        data = sqrt(data.*conj(data));
        temp(jj,:) = data;
    end
    transformBaseline(ii) = mean(temp(:));
end

transformResponse = cell(numChans,2);
for ii=1:numChans
    for jj=1:2
        transformResponse{ii,jj} = zeros(2*reps,stimLen);
        for kk=1:2*reps
            data = Response{ii,jj}(kk,:);
            data = conv(data,kernel,'same');
            data = sqrt(data.*conj(data));
            transformResponse{ii,jj}(kk,:) = data-transformBaseline(ii);
        end
        y = transformResponse{ii,jj};y=y(:);
        if jj==1
            position = horzPosition;
            position = repmat(position,[2*reps,1]);
            
        elseif jj==2
            position = vertPosition;
            position = repmat(position,[2*reps,1]);
        end
        design = [position(:),position(:).*position(:)];
        [b,dev,stats] = glmfit(design,y,'normal','link','log');

        [yhat,dylo,dyhi] = glmval(b,[position(1,:)',position(1,:)'.*position(1,:)'],'log',stats);
        figure();plot(position(:),y,'.');hold on;
        yhat = real(yhat);dylo = real(dylo);dyhi = real(dyhi);
        boundedline(position(1,:)',yhat,[dylo,dyhi],'c');
        title(sprintf('Chan: %d, Dir: %d',ii,jj));
    end
end

% timebandwidth = 60; % approximate standard deviation in time is 
%                 % 0.5*sqrt(timebandwidth/2)
% for ii=1:numChans
%     for jj=1:2
%         figure();
%         for kk=1:2*reps
%             [wt,f] = cwt(Response{ii,jj}(kk,:),sampleFreq,'TimeBandwidth',timebandwidth);
%             wt = sqrt(wt.*conj(wt));
%             [~,ind] = min(abs(f-stimulationFrequency));
%             wt = wt(ind:ind,:);
%             wt = mean(wt,1);
%             baseline = quantile(wt,0.05);
%             subplot(2*reps,1,kk);plot(time,wt-baseline);
%         end
%     end
% end

% discrete or continuous wavelet transform, hilbert transform
%  hilbert(x), dwt, or cwt
% windowLen = floor(checkRefresh*sampleFreq);
% x = 0:windowLen;y = sin(2*pi*x/windowLen);
% z = Response{1,1}(1,:);
% w = conv(z,y);
% plot(w.*w);
end

function [b,dev,se] = GetMLest(design,y,b)




end
