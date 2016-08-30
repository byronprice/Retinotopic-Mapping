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
% Updated: 2016/08/29
%  By: Byron Price

cd('~/CloudStation/ByronExp/Retino');

fileStart = 'RetinoData*.mat';

fileList = dir(fileStart);
numFiles = size(fileList,1);

index = regexp(fileStart,'*');
fileStart = fileStart(1:index-1);
AnimalNames = zeros(numFiles,1);
for ii=1:numFiles
    index = regexp(fileList(ii).name,'_');
    AnimalNames(ii) = str2double(fileList(ii).name(index+1:index+5));
end
AnimalNames = unique(AnimalNames);
numAnimals = length(AnimalNames);

data = struct('Date',cell(numAnimals,1),'Response',cell(numAnimals,1));
datelen = 8;
for ii=1:numAnimals
    list = dir(strcat(fileStart,'*',num2str(AnimalNames(ii)),'*.plx'));
    numDates = length(list);
    data(ii).Dates = zeros(numDates,1);
    data.Images{ii} = cell(numDates,1);
    data.Stats{ii} = cell(numDates,1);
    for jj=1:numDates
        index = regexp(fileList(ii).name,'_');
        Date = str2double(fileList(ii).name(index-datelen:index-1));
        data(ii).Dates(jj,1) = Date;
        [Response] = ExtractSignal(AnimalNames(ii),Date,0);
        % think of way to store Response in a hierarchical structure. It's
        % being output by ExtractSignal as a structure already, so a
        % structure within a structure is really annoying. MATLAB sucks.
    end    
end

% look at the distrubtion of VEP magnitudes on individual trials ... only 
%  in those locations that were deemed significant using MapRetinotopy.m
%  The distribution might be uniform, such that any given VEP magnitude
%  occurs, but I suspect it might be bimodal or at least something close.
%  Fix a Gaussian mixture model and set a threshold between the two the
%  groups. Any VEP magnitude below the threshold will signify, no VEP occurred
%  and anything above will signify a VEP occurred (0 and 1). Then, split up
%  the data in 80% and 20% and train a perceptron on the 80% to take as 
%  input 200 ms of LFP data prior to stimulus onset and try to predict
%  whether or not a VEP will occur. If we've already done the thresholding,
%  then we have the ground truth. Test the model on the 20% of withheld
%  data.

% The question is: is there any information in the LFP that might help us
% determine whether or not a VEP will appear? By setting this up as a
% binary, yes or no, we're effectively getting the perceptron more
% information to work with. If we tried to predict VEP magnitude, then we
% might have only one example of a VEP of size 200 microVolts.

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
    stimAfter = round(0.25*sampleFreq); %250 milliseconds
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
