function [ ] = MappingEffects_Analysis(Animals)
%MappingEffects_Analysis.m
%  Analyze data from mapping effects experiment, see MappingEffects.m
%  Three experimental conditions:
%   1) retinotopic mapping, followed by fake SRP (grey screen)
%   2) retinotopic mapping, followed by SRP
%   3) fake mappin (grey screen), followed by SRP
%
% Created: 2017/02/06
%  Byron Price
% Updated: 2017/02/06
% By: Byron Price

cd('~/CloudStation/ByronExp/MappingEffects');

for ii=1:length(Animals)
   dataFiles = dir(sprintf('MappingEffectsData*%d.plx',Animals(ii)));
   stimFiles = dir(sprintf('MappingEffectsStim*%d.mat',Animals(ii)));
   numFiles = length(dataFiles);numChans = 2;numParameters = 7;stimLen = 300;
   
   
   dailyParameters = zeros(numFiles,numChans,numParameters);
   parameterCI = zeros(numFiles,numChans,numParameters);
   srpVEP = zeros(numFiles,numChans,300);
   
   h = figure();
   for jj=1:numFiles
       ePhysFileName = dataFiles(jj).name;
       [ChanData,timeStamps,tsevs,svStrobed,numChans,sampleFreq] = ExtractSignal(ePhysFileName);
       load(stimFiles(jj).name);
       
       stimLen = round(0.3*sampleFreq); % 300 milliseconds
       vepPosition = round(0.05*sampleFreq):round(0.12*sampleFreq);
       strobeStart = 33;
       strobeTimes = tsevs{strobeStart};
       smoothKernel = 4;
       
       if ConditionNumber == 1
           centerVals = stimParams.centerVals;
           Radius = stimParams.Radius;
           reps = stimParams.reps;
           w_pixels = stimParams.w_pixels;
           h_pixels = stimParams.h_pixels;
           numStimuli = stimParams.numStimuli;
           
           % COLLECT DATA IN THE PRESENCE OF VISUAL STIMULI
           Response = zeros(numChans,numStimuli,reps,stimLen);
           meanResponse = zeros(numChans,numStimuli,stimLen);
           for kk=1:numChans
               for ll=1:numStimuli
                   stimStrobes = strobeTimes(svStrobed == ll);
                   for mm=1:reps
                       stimOnset = stimStrobes(mm);
                       [~,index] = min(abs(timeStamps-stimOnset));
                       temp = ChanData(index:index+stimLen-1,kk);
                       Response(kk,ll,mm,:) = temp;
                   end
                   meanResponse(kk,ll,:) = smooth(mean(squeeze(Response(kk,ll,:,:)),1),smoothKernel);
               end
           end
           
           xaxis = 1:w_pixels;
           yaxis = 1:h_pixels;
           Data = zeros(numChans,numStimuli*reps,2);
           for kk=1:numChans
               count = 1;
               for ll=1:numStimuli
                   for mm=1:reps
                       Data(kk,count,1) = ll;
                       Data(kk,count,2) = min(squeeze(Response(kk,ll,mm,vepPosition)));
                       count = count+1;
                   end
               end
           end
           [finalParameters,fisherInfo,ninetyfiveErrors] = FitLFPGaussRetinoModel(Data,xaxis,yaxis,centerVals);
           MakePlots(finalParameters,meanResponse,xaxis,yaxis,stimLen,Radius,centerVals,numStimuli,numChans,jj,numFiles);
           dailyParameters(jj,:,:) = finalParameters;
           parameterCI(jj,:,:) = ninetyfiveErrors;
       elseif ConditionNumber == 2
           centerVals = stimParams.centerVals;
           Radius = stimParams.Radius;
           reps = stimParams.reps;
           w_pixels = stimParams.w_pixels;
           h_pixels = stimParams.h_pixels;
           numStimuli = stimParams.numStimuli;
           
           % COLLECT DATA IN THE PRESENCE OF VISUAL STIMULI
           Response = zeros(numChans,numStimuli,reps,stimLen);
           meanResponse = zeros(numChans,numStimuli,stimLen);
           for kk=1:numChans
               for ll=1:numStimuli
                   stimStrobes = strobeTimes(svStrobed == ll);
                   for mm=1:reps
                       stimOnset = stimStrobes(mm);
                       [~,index] = min(abs(timeStamps-stimOnset));
                       temp = ChanData(index:index+stimLen-1,kk);
                       Response(kk,ll,mm,:) = temp;
                   end
                   meanResponse(kk,ll,:) = smooth(mean(squeeze(Response(kk,ll,:,:)),1),smoothKernel);
               end
           end
           
           xaxis = 1:w_pixels;
           yaxis = 1:h_pixels;
           Data = zeros(numChans,numStimuli*reps,2);
           for kk=1:numChans
               count = 1;
               for ll=1:numStimuli
                   for mm=1:reps
                       Data(kk,count,1) = ll;
                       Data(kk,count,2) = min(squeeze(Response(kk,ll,mm,vepPosition)));
                       count = count+1;
                   end
               end
           end
           [finalParameters,fisherInfo,ninetyfiveErrors] = FitLFPGaussRetinoModel(Data,xaxis,yaxis,centerVals);
           MakePlots(finalParameters,meanResponse,xaxis,yaxis,stimLen,Radius,centerVals,numStimuli,numChans,jj,numFiles);
           dailyParameters(jj,:,:) = finalParameters;
           parameterCI(jj,:,:) = ninetyfiveErrors;
           
           srpResponse = zeros(numChans,srp_reps,stimLen);
           meanSRP = zeros(numChans,stimLen);
           for kk=1:numChans
               stimStrobes = strobeTimes(svStrobed == srp_word);
               for ll=1:srp_reps
                   stimOnset = stimStrobes(ll);
                   [~,index] = min(abs(timeStamps-stimOnset));
                   temp = ChanData(index:index+stimLen-1,kk);
                   srpResponse(kk,ll,:) = temp;
               end
               meanSRP(kk,:) = smooth(mean(squeeze(srpResponse(kk,:,:)),1),smoothKernel);
               subplot(numFiles,4,ii+2+(jj-1)*4);plot(meanSRP(kk,:));axis([0 stimLen -400 200]);
               xlabel('Time from Phase Reversal (milliseconds)');ylabel('LFP Magnitude (\muVolts)');
               title(sprintf('SRP VEP: Day %d',jj));
           end
           srpVEP(jj,:,:) = meanSRP;
       elseif ConditionNumber == 3
           srpResponse = zeros(numChans,srp_reps,stimLen);
           meanSRP = zeros(numChans,stimLen);
           for kk=1:numChans
               stimStrobes = strobeTimes(svStrobed == srp_word);
               for ll=1:srp_reps
                   stimOnset = stimStrobes(ll);
                   [~,index] = min(abs(timeStamps-stimOnset));
                   temp = ChanData(index:index+stimLen-1,kk);
                   srpResponse(kk,ll,:) = temp;
               end
               meanSRP(kk,:) = smooth(mean(squeeze(srpResponse(kk,:,:)),1),smoothKernel);
               subplot(numFiles,4,kk+2+(jj-1)*4);plot(meanSRP(kk,:));axis([0 stimLen -400 200]);
               if kk==1
                    xlabel('Time from Phase Reversal (milliseconds)');ylabel('LFP Mag (\muVolts)');
                    title(sprintf('SRP VEP: Day %d',jj));
               end
           end
           srpVEP(jj,:,:) = meanSRP;
       end
   end
   savefig(h,sprintf('MappingEffectsResults_%d.fig',Animals(ii)));
   filename = sprintf('MappingEffectsResults_%d.mat',Animals(ii));
   save(filename,'dailyParameters','parameterCI','srpVEP');
end


end

function [ChanData,timeStamps,tsevs,svStrobed,numChans,sampleFreq] = ExtractSignal(EphysFileName)
    % Extract LFP signals from allad, filter, get timestamps
    % read in the .plx file

    if exist(strcat(EphysFileName(1:end-4),'.mat'),'file') ~= 2
        readall(EphysFileName);pause(0.5);
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
%           ChanData(:,ii) = voltage;
    end

    timeStamps = 0:1/sampleFreq:dataLength/sampleFreq-1/sampleFreq;

    if length(timeStamps) ~= dataLength
        display('Error: Review allad cell array and timing')
        return;
    end
    
end

function [] = MakePlots(finalParameters,meanResponse,xaxis,yaxis,stimLen,Radius,centerVals,numStimuli,numChans,Day,numFiles)

    Radius = round(Radius);
    xPos = unique(centerVals(:,1));
    yPos = unique(centerVals(:,2));
    
    xDiff = mean(diff(xPos));
    yDiff = mean(diff(yPos));
    xconv = stimLen/xDiff; 
    yconv = 1000/yDiff; % for height of the stimulus
    
    
    for ii=1:numChans
        subplot(numFiles,4,ii+(Day-1)*4);axis([0 max(xaxis) 0 max(yaxis)]);
        title(sprintf('VEP Retinotopy, Chan %d Day %d',ii,Day));
        xlabel('Horizontal Screen Position (pixels)');ylabel('Vertical Screen Position (pixels)');
        hold on;

        for jj=1:numStimuli
            tempx = centerVals(jj,1);
            tempy = centerVals(jj,2);
            plot(((1:1:stimLen)./xconv+centerVals(jj,1)-0.5*xDiff),...
                (squeeze(meanResponse(ii,jj,:))'./yconv+centerVals(jj,2)),'k','LineWidth',2);

            plot((tempx-Radius)*ones(Radius*2+1,1),(tempy-Radius):(tempy+Radius),'k','LineWidth',2);
            plot((tempx+Radius)*ones(Radius*2+1,1),(tempy-Radius):(tempy+Radius),'k','LineWidth',2);
            plot((tempx-Radius):(tempx+Radius),(tempy-Radius)*ones(Radius*2+1,1),'k','LineWidth',2);
            plot((tempx-Radius):(tempx+Radius),(tempy+Radius)*ones(Radius*2+1,1),'k','LineWidth',2);

        end
        
        finalIm = zeros(length(xaxis),length(yaxis));
        parameterVec = finalParameters(ii,:);
        for jj=1:length(x)
            for kk=1:length(y)
                distX = x(jj)-parameterVec(2);
                distY = y(kk)-parameterVec(3);
                b = [parameterVec(1),parameterVec(4),parameterVec(5),parameterVec(7)];
                finalIm(jj,kk) = b(1)*exp(-(distX.^2)./(2*b(2)*b(2))-(distY.^2)./(2*b(3)*b(3)))+b(4);
            end
        end
        imagesc(x,y,finalIm','AlphaData',0.5);set(gca,'YDir','normal');w=colorbar;
        ylabel(w,'Mean Single-Trial VEP Negativity (\muV)');colormap('jet');hold off;
    end
end


