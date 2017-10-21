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
%       Chans - channel numbers, input as [6,8], defaults to 6 and 8
%OUTPUT: plots
%
% Created: 2016/05/31, 24 Cummington, Boston
%  Byron Price
% Updated: 2016/08/22
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
Response = cell(numChans,numDirs);

for ii=1:numChans
    for jj=1:numDirs
        strobeNum = jj;
        currentStrobeTimes = strobeTimes(svStrobed==strobeNum);
        Response{ii,jj} = zeros(reps,stimLen);
        for kk=1:reps
           onsetTime = round(currentStrobeTimes(kk)*sampleFreq);
           offsetTime = onsetTime+stimLen-1;
           Response{ii,jj}(kk,:) = ChanData(onsetTime:offsetTime,ii);
        end
    end
end

% create sliding window spectrogram for power at the stimulation frequency as function
%  of position on the screen
fftSpectrogram = cell(numChans,numDirs);

windowLen = floor(0.5*sampleFreq);
nOverlap = floor(windowLen/5);
nfft = max(256,2^nextpow2(windowLen));
for ii=1:numChans
    for jj=1:numDirs
        powerAtStimFreq = zeros(reps,5,
        for kk=1:reps
            data = Response{ii,jj}(kk,:)';
            [s,f,t] = spectrogram(data,windowLen,nOverlap,nfft,sampleFreq);
            [~,ind] = min(abs(f-stimulationFrequency));
            tempfft = s(ind-2:ind+2,:);
            power = sqrt(tempfft.*conj(tempfft));
        end
    end
end
pause(1);
end
