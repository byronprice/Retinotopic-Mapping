function [] = MapRetinotopy(AnimalName,Date)
% MapRetinotopy.m
%
%  Will take data from a retinotopic mapping experiment and extract the
%   retinotopy of the LFP recording electrode.
%INPUT: AnimalName - unique identifier for the animal as a number, e.g.
%            12345
%       Date - date of the experiment, e.g. 20160525

%OUTPUT: saved files with info regarding retinotopy of each
%         channel ... parameters defined elsewhere in
%         FitLFPRetinoModel_LM.m
%
% Created: 2016/05/25, 8 St. Mary's Street, Boston
%  Byron Price
% Updated: 2017/02/21
%  By: Byron Price

cd('~/CloudStation/ByronExp/Retino');

EphysFileName = sprintf('RetinoData%d_%d.plx',Date,AnimalName); 
  
global centerVals Radius stimTime holdTime numStimuli w_pixels h_pixels ...
    DistToScreen numChans sampleFreq stimLen minWin maxWin stimStrobeNum; %#ok<*REDEF>

StimulusFileName = sprintf('RetinoStim%d_%d.mat',Date,AnimalName);
load(StimulusFileName)

fprintf('Opening File: %s ...\n',EphysFileName);

centerVals = stimParams.centerVals;
Radius = stimParams.Radius;
stimStrobeNum = stimParams.stimStrobeNum;
stimTime = stimParams.stimTime;
holdTime = stimParams.holdTime;
numStimuli = stimParams.numStimuli;
w_pixels = stimParams.w_pixels;
h_pixels = stimParams.h_pixels;
DistToScreen = stimParams.DistToScreen;
mmPerPixel = stimParams.mmPerPixel;

% convert from allad to ChanData by filtering
[ChanData,timeStamps,tsevs,svStrobed] = ExtractSignal(EphysFileName);

% get LFP response to each stimulus (the VEPs)
stimLen = round(0.3*sampleFreq);
minWin = round(0.05*sampleFreq):1:round(0.12*sampleFreq);
maxWin = round(.15*sampleFreq):1:round(0.3*sampleFreq);
[Response] = CollectVEPS(ChanData,timeStamps,tsevs,svStrobed);

xaxis = 1:w_pixels;
yaxis = 1:h_pixels;
Data = cell(numChans,1);

for ii=1:numChans
    tempData = zeros(numStimuli,2);
    for jj=1:numStimuli
        tempData(jj,1) = jj;
        tempData(jj,2) = max(squeeze(Response(ii,jj,maxWin)))-min(squeeze(Response(ii,jj,minWin)));
    end
    temp = tempData(:,2);
    outlier = mean(temp)+4*std(temp);
    indeces = find(temp>outlier);
    tempData(indeces,:) = [];
    Data{ii} = tempData;
end

fprintf('Fitting model ...\n');
[finalParameters,fisherInfo,ninetyfiveErrors] = FitLFPRetinoModel_LM(Data,xaxis,yaxis,centerVals);

fprintf('Making plots ...\n');
[h] = MakePlots(finalParameters,AnimalName,xaxis,yaxis); 

vepResponse = Response;
dimReduceData = Data;
save(sprintf('RetinoMap%d_%d.mat',Date,AnimalName),'vepResponse','dimReduceData','finalParameters','fisherInfo','ninetyfiveErrors',...
    'numChans','w_pixels','h_pixels','mmPerPixel');
end

function [ChanData,timeStamps,tsevs,svStrobed] = ExtractSignal(EphysFileName)
    % Extract LFP signals from allad, filter, get timestamps
    global numChans sampleFreq;
    % read in the .plx file

    if exist(strcat(EphysFileName(1:end-4),'.mat'),'file') ~= 2
        readall(EphysFileName);
    end

    EphysFileName = strcat(EphysFileName(1:end-4),'.mat');
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
        fprintf('Error: Review allad cell array and timing')
        return;
    end
    
end

function [Response] = CollectVEPS(ChanData,timeStamps,tsevs,svStrobed)
    global numChans numStimuli stimLen stimStrobeNum;
    strobeStart = 33;
    strobeTimes = tsevs{1,strobeStart};
    % COLLECT DATA IN THE PRESENCE OF VISUAL STIMULI

    Response = zeros(numChans,numStimuli,stimLen);

    for ii=1:numChans
        stimStrobes = strobeTimes(svStrobed == stimStrobeNum);
        if numStimuli == length(stimStrobes)
            for jj=1:numStimuli
                stimOnset = stimStrobes(jj);
                [~,index] = min(abs(timeStamps-stimOnset));
                temp = ChanData(index:index+stimLen-1,ii);
                Response(ii,jj,:) = temp;
            end
        else 
            fprintf('Error: Number of strobes does not equal number of stimuli');
        end
    end
end

function [h] = MakePlots(finalParameters,AnimalName,x,y)
    global numChans w_pixels h_pixels;

    
    for ii=1:numChans
        h(ii) = figure;
    end

    for ii=1:numChans
        figure(h(ii));axis([0 w_pixels 0 h_pixels]);
        title(sprintf('VEP Retinotopy: Chan %d, Animal %d',ii,AnimalName));
        xlabel('Horizontal Screen Position (pixels)');ylabel('Vertical Screen Position (pixels)');
        hold on;
        
        finalIm = zeros(length(x),length(y));
        parameterVec = finalParameters(ii,:);
        for jj=1:length(x)
            for kk=1:length(y)
                distX = x(jj)-parameterVec(2);
                distY = y(kk)-parameterVec(3);
                b = [parameterVec(1),parameterVec(4),parameterVec(5),parameterVec(7)];
                finalIm(jj,kk) = b(1)*exp(-(distX.^2)./(2*b(2)*b(2))-(distY.^2)./(2*b(3)*b(3)))+b(4);
            end
        end
        imagesc(x,y,finalIm');set(gca,'YDir','normal');w=colorbar;
        caxis([parameterVec(7) parameterVec(7)+500]);
        ylabel(w,'Mean VEP Magnitude (\muV)');colormap('jet');hold off;
        
    end
end
