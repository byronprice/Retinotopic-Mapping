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
% Updated: 2016/06/08
%  By: Byron Price

EphysFileName = strcat('RetinoData',num2str(Date),'_',num2str(AnimalName));

readall(EphysFileName);

StimulusFileName = strcat('RetinoStim',num2str(Date),'_',num2str(AnimalName),'.mat');
EphysFileName = strcat(EphysFileName,'.mat');
load(EphysFileName)
load(StimulusFileName)

if nargin < 3
    Chans = [6,8];
end

sampleFreq = adfreq;

% tsevs are the strobed times of stimulus onset, then offset
%  Onset at tsevs{1,33}(2), offset at tsevs{1,33}(3), onset at
%  tsevs{1,33}(4), offset at 5, etc.
% allad contains the continuous data from each channel, which appear to be
%  recorded at 1000 Hz rather than 40,000
x = find(~cellfun(@isempty,tsevs));
strobeStart = x(1);

dataLength = length(allad{1,strobeStart+Chans(1)-1});
numChans = length(Chans);
ChanData = zeros(dataLength,numChans);
for ii=1:numChans
    ChanData(:,ii) = allad{1,strobeStart+Chans(ii)-1};
end
timeStamps = 0:1/sampleFreq:dataLength/sampleFreq-1/sampleFreq;

if length(timeStamps) ~= dataLength
    display('Error: Review allad cell array and timing')
    return;
end
strobeData = tsevs{1,strobeStart};
totalStrobes = length(strobeData);

if mod(totalStrobes,2) == 1
    display('Error: Missed a stimulus onset/offset strobe')
    return;
end

reps = 5;
dataPoints = length(strobeData(1:2:end-1));
B = glmfit(1:dataPoints,strobeData(2:2:end)-strobeData(1:2:end-1),'Normal');
timeWindow = B(1);

Response = zeros(dataPoints/reps,numChans);
timeFrames = round(timeWindow*sampleFreq);

delay = round(0.04*sampleFreq); % add delay after stimulus onset to account for
         % ~40ms delay in neuronal processing
for ii=1:numChans
    count = 1;
    for jj=1:dataPoints/reps
        temp = 0;
        for kk=1:reps
            stimOnset = strobeData(jj+dataPoints/reps*(kk-1));
            [~,index] = min(abs(timeStamps-stimOnset));
            temp = temp+ChanData(index+delay:index+delay+timeFrames,ii);
        end
        avg = temp./reps;
        Response(count,ii) = max(avg)-min(avg);
        count = count+1;
    end
end
cutOff = prctile(Response,98);

Indeces = cell(1,numChans);
figure();hold on;
for ii=1:numChans
        Indeces{ii} = find(Response(:,ii)>cutOff(ii));
        subplot(2,1,ii);plot(stimulusLocs(Indeces{ii},1),stimulusLocs(Indeces{ii},2),'*b','LineWidth',2);
end
hold off;

end
