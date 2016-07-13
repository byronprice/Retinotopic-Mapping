function [] = MapRetinotopy(AnimalName,Date,Chans)
% MapRetinotopy.m
%
%  Will take data from a retinotopic mapping experiment and extract the
%   retinotopy of the LFP recording electrode.
%INPUT: AnimalName - unique identifier for the animal as a number, e.g.
%            12345
%       Date - date of the experiment, e.g. 20160525
%       Chans - channel numbers, input as [6,8], defaults to 6 and 8
%OUTPUT:
%
% Created: 2016/05/25, 8 St. Mary's Street, Boston
%  Byron Price
% Updated: 2016/07/13
%  By: Byron Price


% read in the .plx file
EphysFileName = strcat('RetinoData',num2str(Date),'_',num2str(AnimalName));

if exist(strcat(EphysFileName,'.mat'),'file') ~= 2
    MyReadall(EphysFileName);
end

StimulusFileName = strcat('RetinoStim',num2str(Date),'_',num2str(AnimalName),'.mat');
EphysFileName = strcat(EphysFileName,'.mat');
load(EphysFileName)
load(StimulusFileName)

if nargin < 3
    Chans = [6,8];
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

%x = find(~cellfun(@isempty,tsevs));
strobeStart = 33;

% lowpass filter the data
dataLength = length(allad{1,strobeStart+Chans(1)-1});
numChans = length(Chans);
ChanData = zeros(dataLength,numChans);
preAmpGain = 1;
for ii=1:numChans
    voltage = ((allad{1,strobeStart+Chans(ii)-1}).*SlowPeakV)./(0.5*(2^SlowADResBits)*adgains(strobeStart+Chans(ii)-1)*preAmpGain);
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
stimLength = round((stimLen+0.5)*sampleFreq/2)*2;

% COLLECT DATA IN THE PRESENCE OF VISUAL STIMULI
Response = zeros(numChans,numStimuli,reps,stimLength);
for ii=1:numChans
    for jj=1:reps
        check = (jj-1)*numStimuli+1:jj*numStimuli;
        for kk=1:numStimuli
            stimOnset = strobeData(check(kk));
            [~,index] = min(abs(timeStamps-stimOnset));
            temp = ChanData(index:index+stimLength-1,ii);
            Response(ii,kk,jj,:) = temp;
        end
        clear check;
    end
end

meanResponse = squeeze(mean(Response,3));
Statistic = abs(min(meanResponse,[],3));

% STATISTIC OF INTEREST is T = abs(min(mean(LFP across stimulus repetitions)) 
% in the interval from 0 to ~ 0.5 seconds after an image is flashed on the 
% screen, this is a measure of the size of a VEP

% BOOTSTRAP FOR 95% CONFIDENCE INTERVALS OF STATISTIC IN ABSENCE OF VISUAL STIMULI
%  started trial with 30 seconds of a blank screen
N = 5000; % number of bootstrap samples
noStimLen = startPause*sampleFreq;

bootPrctile = zeros(numChans,1); % 99 percentile
for ii=1:numChans
    Tboot = zeros(N,1);
    for jj=1:N
        indeces = random('Discrete Uniform',noStimLen-stimLength*2,[reps,1]);
        temp = zeros(reps,stimLength);
        for kk=1:reps
            temp(kk,:) = ChanData(indeces(kk):indeces(kk)+stimLength-1,ii);
        end
        Tboot(jj) = abs(min(mean(temp,1)));
    end
    %figure();histogram(Tboot);
    bootPrctile(ii) = quantile(Tboot,1-1/100);
end

for ii=1:numChans
    figure();histogram(Statistic(ii,:));
    hold on; plot(bootPrctile(ii)*ones(100,1),0:99,'LineWidth',2);
end

significantStimuli = zeros(numChans,numStimuli);
for ii=1:numChans
    for jj=1:numStimuli
        if Statistic(ii,jj) >= bootPrctile(ii)
            significantStimuli(ii,jj) = Statistic(ii,jj)./bootPrctile(ii);
        end
    end    
end

stimVals = zeros(numChans,w_pixels,h_pixels);
x=1:w_pixels;
y=1:h_pixels;
for ii=1:numChans
    for jj=1:numStimuli
        tempx = centerVals(jj,1);
        tempy = centerVals(jj,2);
%         for kk=(tempx-Radius):(tempx+Radius)
%             for ll=(tempy-Radius):(tempy+Radius)
%                 pointx = kk-tempx;
%                 pointy = ll-tempy;
%                 if pointx*pointx+pointy*pointy <= Radius*Radius
%                     stimVals(ii,kk,ll) = significantStimuli(ii,jj);
%                 end
%             end
%         end
        stimVals(ii,tempx-Radius:tempx+Radius,tempy-Radius:tempy+Radius) = significantStimuli(ii,jj);
    end
    figure();imagesc(x,y,squeeze(stimVals(ii,:,:))');set(gca,'YDir','normal');colorbar;
end

end
