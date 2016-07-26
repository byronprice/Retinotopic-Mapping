function [stimVals,centerMass,numChans] = MapRetinotopy(AnimalName,Date,yesNo)
% MapRetinotopy.m
%
%  Will take data from a retinotopic mapping experiment and extract the
%   retinotopy of the LFP recording electrode.
%INPUT: AnimalName - unique identifier for the animal as a number, e.g.
%            12345
%       Date - date of the experiment, e.g. 20160525
%       yesNo - 1 if you are running this code as is, 0 if it is being run
%        through MapRetWrapper.m
%OUTPUT: stimVals
%        x -
%        y -
%        centerMass - 
%
% Created: 2016/05/25, 8 St. Mary's Street, Boston
%  Byron Price
% Updated: 2016/07/25
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
preAmpGain = 1;
for ii=1:numChans
    voltage = 1000.*((allad{1,Chans(ii)}).*SlowPeakV)./(0.5*(2^SlowADResBits)*adgains(Chans(ii))*preAmpGain);
    temp = smooth(voltage,0.013*sampleFreq);
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
strobeTimes = tsevs{1,strobeStart};
stimLength = round((stimLen+0.2)*sampleFreq); % about 250 milliseconds

% COLLECT DATA IN THE PRESENCE OF VISUAL STIMULI
Response = zeros(numChans,numStimuli,reps,stimLength);
for ii=1:numChans
    for jj=1:reps
        check = (jj-1)*numStimuli+1:jj*numStimuli;
        for kk=1:numStimuli
            stimOnset = strobeTimes(check(kk));
            [~,index] = min(abs(timeStamps-stimOnset));
            temp = ChanData(index:index+stimLength-1,ii);
            Response(ii,kk,jj,:) = temp;
        end
        clear check;
    end
end

% STATISTIC OF INTEREST is T = max - min(mean(LFP across stimulus repetitions)) 
% in the interval from 0 to ~ 0.3 seconds after an image is flashed on the 
% screen, this is a measure of the size of a VEP
meanResponse = squeeze(mean(Response,3));
dataStat = max(meanResponse,[],3)-min(meanResponse,[],3);

% BOOTSTRAP FOR STANDARD ERROR OF STATISTIC IN PRESENCE OF VISUAL STIMULI
N = 5000;
dataError = zeros(numChans,numStimuli);
for ii=1:numChans
    for jj=1:numStimuli
        Tboot = zeros(N,1);
        for kk=1:N
            indeces = random('Discrete Uniform',reps,[reps,1]);
            group = squeeze(Response(ii,jj,indeces,:));
            meanGroup = mean(group,1);
            Tboot(kk) = max(meanGroup)-min(meanGroup);
        end
        dataError(ii,jj) = std(Tboot);
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
        indeces = random('Discrete Uniform',noStimLen-stimLength*2,[reps,1]);
        temp = zeros(reps,stimLength);
        for kk=1:reps
            temp(kk,:) = ChanData(indeces(kk):indeces(kk)+stimLength-1,ii);
        end
        Tboot(jj) = max(mean(temp,1))-min(mean(temp,1));
    end
    %figure();histogram(Tboot);
    bootPrctile(ii) = quantile(Tboot,1-1/100);
    bootStat(ii,1) = mean(Tboot);
    bootStat(ii,2) = std(Tboot);
end

% for ii=1:numChans
%     figure();histogram(dataStat(ii,:));
%     hold on; plot(bootPrctile(ii)*ones(100,1),0:99,'LineWidth',2);
% end

% WALD TEST - VEP magnitude significantly greater in presence of a stimulus
%  than in the absence of a stimulus
significantStimuli = zeros(numChans,numStimuli);
alpha = 0.05/numStimuli;
for ii=1:numChans
    for jj=1:numStimuli
        W = (dataStat(ii,jj)-bootStat(ii,1))/sqrt(dataError(ii,jj)^2+bootStat(ii,2)^2);
        c = norminv(1-alpha,0,1);
        if W > c
            significantStimuli(ii,jj) = dataStat(ii,jj); % or equals W itself
        end
    end    
end

stimVals = zeros(numChans,w_pixels,h_pixels);
x=1:w_pixels;
y=1:h_pixels;
if yesNo == 1
    figure();
end
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
    if yesNo == 1
        subplot(numChans,1,ii);imagesc(x,y,squeeze(stimVals(ii,:,:))');set(gca,'YDir','normal');h=colorbar;
        title(sprintf('Retinotopy for Channel %d , Animal %d',ii,AnimalName));ylabel(h,'VEP Magnitude (\muV)');
        xlabel('Horizontal Screen Position (pixels)');ylabel('Vertical Screen Position (pixels)');
        colormap('jet');
    end
end
if yesNo == 1
    savefig(strcat('RetinoMap',num2str(Date),'_',num2str(AnimalName),'.fig'));
end

centerMass = zeros(numChans,4);
Sigma = zeros(numChans,2,2);
for ii=1:numChans
    dataX = [];dataY = [];
    for jj=1:numStimuli
        dataX = [dataX;repmat(centerVals(jj,1),[round(significantStimuli(ii,jj)),1])];
        dataY = [dataY;repmat(centerVals(jj,2),[round(significantStimuli(ii,jj)),1])];
    end
    data = [dataX,dataY];

    try
        mnPDF = fitgmdist(data,1);
        centerMass(ii,1) = mnPDF.mu(1);
        centerMass(ii,2) = mnPDF.mu(2);
        centerMass(ii,3) = mnPDF.Sigma(1,1);
        centerMass(ii,4) = mnPDF.Sigma(2,2);
        Sigma(ii,:,:) = mnPDF.Sigma;
    catch exception
        display(sprintf('Error on Channel %d. Do not trust centerMass values.',ii));
        centerMass(ii,:) = NaN;
        Sigma(ii,:,:) = NaN;
    end

end

save(strcat('RetinoMap',num2str(Date),'_',num2str(AnimalName),'.mat'),'numChans',...
    'centerVals','significantStimuli','centerMass','stimVals','Sigma');

% obj = gmdistribution(centerMass(chan,1:2),squeeze(Sigma(chan,:,:)));
% figure();
% h = ezcontour(@(x,y) pdf(obj,[x y]),[0 w_pixels,[0 h_pixels]);
end
