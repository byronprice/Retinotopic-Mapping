function [] = PredictVEP()
% PredictVEP.m
%
%  Will take data from all retinotopy exps and run a machine learning 
%   algorithm to see if you can predict when a VEP will occur
%INPUT: 
%OUTPUT: 
%
% Created: 2016/08/02
%  Byron Price
% Updated: 2016/08/02
%  By: Byron Price

cd('~/CloudStation/ByronExp/RetinoExp');

fileStart = 'RetinoData*.mat';

fileList = dir(fileStart);
numFiles = size(fileList,1);

fileStart = fileList(1).name(1:end-18);
AnimalNames = zeros(numFiles,1);
for ii=1:numFiles
    AnimalNames(ii) = str2double(fileList(ii).name(end-8:end-4));
end
AnimalNames = unique(AnimalNames);
numAnimals = length(AnimalNames);

data = struct('Date',cell(numAnimals,1),'Response',cell(numAnimals,1));
for ii=1:numAnimals
    list = dir(strcat(fileStart,'*',num2str(AnimalNames(ii)),'*.plx'));
    numDates = length(list);
    data(ii).Dates = zeros(numDates,1);
    data.Images{ii} = cell(numDates,1);
    data.Stats{ii} = cell(numDates,1);
    for jj=1:numDates
        Date = str2double(list(jj).name(end-17:end-10));
        data(ii).Dates(jj,1) = Date;
        [Response] = ExtractSignal(AnimalNames(ii),Date,0);
        
    end    
end


end

function [Response] = ExtractSignal(AnimalName,Date)
% Get LFP in window before and after each stimulus onset

    % read in the .plx file
    EphysFileName = strcat('RetinoData',num2str(Date),'_',num2str(AnimalName));

    if exist(strcat(EphysFileName,'.mat'),'file') ~= 2
        MyReadall(EphysFileName);
    end

    StimulusFileName = strcat('RetinoStim',num2str(Date),'_',num2str(AnimalName),'.mat');
    EphysFileName = strcat(EphysFileName,'.mat');
    load(EphysFileName)
    load(StimulusFileName)

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
    stimAfter = round((stimTime+0.2)*sampleFreq); % about 250 milliseconds
    stimBefore = round(0.1*sampleFreq); %100 milliseconds

    % COLLECT DATA IN THE PRESENCE OF VISUAL STIMULI
    Response = struct('Before',zeros(numChans,numStimuli,reps,stimBefore),...
        'After',zeros(numChans,numStimuli,reps,stimAfter));
    for ii=1:numChans
        for jj=1:numStimuli
            stimStrobes = strobeTimes(svStrobed == jj);
            for kk=1:reps
                stimOnset = stimStrobes(kk);
                [~,index] = min(abs(timeStamps-stimOnset));
                temp1 = ChanData(index-stimBefore+1:index);
                temp2 = ChanData(index:index+stimAfter-1,ii);
                Response.Before(ii,jj,kk,:) = temp1;
                Response.After(ii,jj,kk,:) = temp2;
            end
            clear temp1 temp2;
        end
    end


end
