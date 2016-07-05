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
% Updated: 2016/07/05
%  By: Byron Price


% read in the .plx file
EphysFileName = strcat('RetinoData',num2str(Date),'_',num2str(AnimalName));

if exist(strcat(EphysFileName,'.mat'),'file') ~= 2
    readall(EphysFileName);
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
for ii=1:numChans
    temp = smooth(allad{1,strobeStart+Chans(ii)-1},0.05*sampleFreq);
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

% BOOTSTRAP FOR 95% CONFIDENCE INTERVALS OF STATISTIC IN ABSENCE OF VISUAL STIMULI
%  started trial with 30 seconds of a blank screen
N = 1000; % number of bootstrap samples
noStimLen = startPause*sampleFreq;

bootPrctile = zeros(numChans,1); % 99.0 percentile
for ii=1:numChans
    Tboot = zeros(N,1);
    indeces = randperm(noStimLen-stimLength*2,N);
    for jj=1:N
        temp = ChanData(indeces(jj):indeces(jj)+stimLength-1,ii);
        Tboot(jj) = max(temp)-min(temp);
    end
    figure();histogram(Tboot);
    bootPrctile(ii) = quantile(Tboot,0.99);
end

% CALCULATE STATISTIC IN PRESENCE OF VISUAL STIMULI
Response = zeros(numChans,numStimuli,reps);
for ii=1:numChans
    for jj=1:reps
        check = (jj-1)*numStimuli+1:jj*numStimuli;
        for kk=1:numStimuli
            stimOnset = strobeData(check(kk));
            [~,index] = min(abs(timeStamps-stimOnset));
            temp = ChanData(index:index+stimLength-1,ii);
            Response(ii,kk,jj) = max(temp)-min(temp);
        end
        clear check;
    end
end

for ii=1:numChans
    figure();histogram(reshape(squeeze(Response(ii,:,:)),[numStimuli*reps,1]));
    hold on; plot(bootPrctile(ii)*ones(100,1),0:99,'LineWidth',2);
end

significantStimuli = zeros(numChans,numStimuli);
for ii=1:numChans
    for jj=1:numStimuli
        for kk=1:reps
            if Response(ii,jj,kk) >= bootPrctile(ii)
                significantStimuli(ii,jj) = significantStimuli(ii,jj)+1./reps;
            end
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
        for kk=(tempx-Radius):(tempx+Radius)
            pointx = (kk-tempx)^2;
            for ll=(tempy-Radius):(tempy+Radius)
                pointy = (ll-tempy)^2;
                if (pointx+pointy) <= (Radius*Radius)
                    stimVals(ii,kk,ll) = significantStimuli(ii,jj);
                end
            end
        end
    end
    figure();imagesc(x,y,squeeze(stimVals(ii,:,:)));set(gca,'YDir','normal');colorbar;
end

end
