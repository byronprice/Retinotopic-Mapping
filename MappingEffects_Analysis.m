function [ ] = MappingEffects_Analysis(Animals,Channels)
%MappingEffects_Analysis.m
%  Analyze data from mapping effects experiment, see MappingEffects.m
%  Three experimental conditions:
%   1) retinotopic mapping, followed by fake SRP (grey screen)
%   2) retinotopic mapping, followed by SRP
%   3) fake mapping (grey screen), followed by SRP
%
% Created: 2017/02/06
%  Byron Price
% Updated: 2017/03/27
% By: Byron Price

cd('~/CloudStation/ByronExp/MappingEffects');

for ii=1:length(Animals)
   dataFiles = dir(sprintf('MappingEffectsData*%d.plx',Animals(ii)));
   stimFiles = dir(sprintf('MappingEffectsStim*%d.mat',Animals(ii)));
   
   numFiles = length(dataFiles);
   fprintf('Running animal %d\n\n',Animals(ii));
   
   pix_to_degree = cell(numFiles,1);
   dailyParameters = cell(numFiles,1);
   mapData = cell(numFiles,1);
   parameterCI = cell(numFiles,1);
   fisherInformation = cell(numFiles,1);
   srpSize = cell(numFiles,1);
   srpVEP = cell(numFiles,1);
   residualDevianceTest_pVal = cell(numFiles,1);
   fullDeviance = cell(numFiles,1);
   significantVEP = cell(numFiles,1);
   currentChannel = Channels{ii};
   numChans = length(currentChannel);
   mappingVEP = cell(numFiles,numChans);
   srpTrialVEPs = cell(numFiles,numChans);
   
   % could add back check for "goodChannels"
   if numChans ~= 0
       h = figure();
       for jj=1:numFiles
           ePhysFileName = dataFiles(jj).name;
           [ChanData,timeStamps,tsevs,svStrobed,~,sampleFreq] = ExtractSignal(ePhysFileName);
           load(stimFiles(jj).name);
           
           stimLen = round(0.5*sampleFreq); % 500 milliseconds
           strobeStart = 33;
           strobeTimes = tsevs{strobeStart};
           
           
           if ConditionNumber == 1
               centerVals = stimParams.centerVals;
               Radius = stimParams.Radius;
               reps = stimParams.reps;
               w_pixels = stimParams.w_pixels;
               h_pixels = stimParams.h_pixels;
               numStimuli = stimParams.numStimuli;
               mmPerPixel = stimParams.mmPerPixel;
               DistToScreen = stimParams.DistToScreen*10;
               
               pix_to_degree{jj} = mmPerPixel/DistToScreen;
               
               
               % COLLECT DATA IN THE PRESENCE OF VISUAL STIMULI
               Response = zeros(numChans,numStimuli,reps,stimLen);
               meanResponse = zeros(numChans,numStimuli,stimLen);
               for kk=1:numChans
                   for ll=1:numStimuli
                       stimStrobes = strobeTimes(svStrobed == ll);
                       for mm=1:reps
                           stimOnset = stimStrobes(mm);
                           [~,index] = min(abs(timeStamps-stimOnset));
                           temp = ChanData(index:index+stimLen-1,currentChannel(kk));
                           Response(kk,ll,mm,:) = temp;
                       end
                       meanResponse(kk,ll,:) = mean(squeeze(Response(kk,ll,:,:)),1);
                   end
               end

               xaxis = 1:w_pixels;
               yaxis = 1:h_pixels;
               Data = cell(numChans,1);
               for kk=1:numChans
                   count = 1;Data{kk} = zeros(reps*numStimuli,3);
                   
                   allVEPs = reshape(squeeze(Response(kk,:,:,:)),[reps*numStimuli,stimLen]);
                   meanVEP = mean(allVEPs,1);
                   mappingVEP{jj,kk} = allVEPs;
                  
                   [~,minLatency] = min(meanVEP);
                
                   minWin = (minLatency-30):(minLatency+30);
                   maxWin = (minLatency+40):(minLatency+140);
                 
                   
                   for ll=1:numStimuli
                       for mm=1:reps
                           Data{kk}(count,1:2) = centerVals(ll,:);
                           [minVal,~] = min(squeeze(Response(kk,ll,mm,minWin)));
                           Data{kk}(count,3) = max(squeeze(Response(kk,ll,mm,maxWin)))...
                               -minVal;
                           count = count+1;
                       end
                   end
                   tempData = Data{kk};
                   temp = tempData(:,3);
                   outlier = median(temp)+10*std(temp);
                   indeces = find(temp>outlier);
                   tempData(indeces,:) = [];
                   
                   temp = tempData(:,3);
                   indeces = find(temp==0);
                   tempData(indeces,:) = [];
                   Data{kk} = abs(tempData);
               end

               numRepeats = 5e3;
               [finalParameters,fisherInfo,ninetyfiveErrors,signifMap,Deviance,residDevTestp] = ...
                        FitLFPRetinoModel_Gamma(Data,xaxis,yaxis,numRepeats);
               MakePlots(finalParameters,meanResponse,xaxis,yaxis,stimLen,Radius,centerVals,numStimuli,numChans,jj,numFiles,h,ConditionNumber);
               
%                x = input('Maps okay? (y/n): ','s');
%                if x == 'n'
%                    numRepeats = 1e4;
%                    [finalParameters,fisherInfo,ninetyfiveErrors,signifMap,Deviance,residDevTestp] = ...
%                         FitLFPRetinoModel_Loglog(Data,xaxis,yaxis,numRepeats);
%                         MakePlots(finalParameters,meanResponse,xaxis,yaxis,stimLen,Radius,centerVals,numStimuli,numChans,jj,numFiles,h,ConditionNumber);
%                end
               
               dailyParameters{jj} = finalParameters;
               mapData{jj} = Data;
               parameterCI{jj} = ninetyfiveErrors;
               residualDevianceTest_pVal{jj} = residDevTestp;
               fullDeviance{jj} = Deviance;
               fisherInformation{jj} = fisherInfo;
           elseif ConditionNumber == 2
               centerVals = stimParams.centerVals;
               Radius = stimParams.Radius;
               reps = stimParams.reps;
               w_pixels = stimParams.w_pixels;
               h_pixels = stimParams.h_pixels;
               numStimuli = stimParams.numStimuli;
               mmPerPixel = stimParams.mmPerPixel;
               DistToScreen = stimParams.DistToScreen*10;
               
               pix_to_degree{jj} = mmPerPixel/DistToScreen;
               
               % COLLECT DATA IN THE PRESENCE OF VISUAL STIMULI
               Response = zeros(numChans,numStimuli,reps,stimLen);
               meanResponse = zeros(numChans,numStimuli,stimLen);
               for kk=1:numChans
                   for ll=1:numStimuli
                       stimStrobes = strobeTimes(svStrobed == ll);
                       for mm=1:reps
                           stimOnset = stimStrobes(mm);
                           [~,index] = min(abs(timeStamps-stimOnset));
                           temp = ChanData(index:index+stimLen-1,currentChannel(kk));
                           Response(kk,ll,mm,:) = temp;
                       end
                       meanResponse(kk,ll,:) = mean(squeeze(Response(kk,ll,:,:)),1);
                   end
               end
               
               xaxis = 1:w_pixels;
               yaxis = 1:h_pixels;
               Data = cell(numChans,1);
               for kk=1:numChans
                   count = 1;Data{kk} = zeros(reps*numStimuli,3);
                   
                   allVEPs = reshape(squeeze(Response(kk,:,:,:)),[reps*numStimuli,stimLen]);
                   meanVEP = mean(allVEPs,1);
                   mappingVEP{jj,kk} = allVEPs;
                   
                   [~,minLatency] = min(meanVEP);

                   minWin = (minLatency-30):(minLatency+30);
                   maxWin = (minLatency+40):(minLatency+140);

                   
                   for ll=1:numStimuli
                       for mm=1:reps
                           Data{kk}(count,1:2) = centerVals(ll,:);
                           [minVal,~] = min(squeeze(Response(kk,ll,mm,minWin)));
                           Data{kk}(count,3) = max(squeeze(Response(kk,ll,mm,maxWin)))...
                               -minVal;
                           count = count+1;
                       end
                   end
                   tempData = Data{kk};
                   temp = tempData(:,3);
                   outlier = median(temp)+10*std(temp);
                   indeces = find(temp>outlier);
                   tempData(indeces,:) = [];
                   
                   temp = tempData(:,3);
                   indeces = find(temp==0);
                   tempData(indeces,:) = [];
                   Data{kk} = abs(tempData);
               end

%                [finalParameters,fisherInfo,ninetyfiveErrors] = FitLFPRetinoModel_Gamma(Data,xaxis,yaxis,centerVals);
%                [finalParameters,covariance] = BayesianFitLFPModel(Data,xaxis,yaxis,centerVals);
               numRepeats = 5e3;
               [finalParameters,fisherInfo,ninetyfiveErrors,signifMap,Deviance,residDevTestp] = ...
                        FitLFPRetinoModel_Gamma(Data,xaxis,yaxis,numRepeats);
                    
               MakePlots(finalParameters,meanResponse,xaxis,yaxis,stimLen,Radius,centerVals,numStimuli,numChans,jj,numFiles,h,ConditionNumber);
               
%                x = input('Maps okay? (y/n): ','s');
%                if x == 'n'
%                    numRepeats = 1e4;
%                    [finalParameters,fisherInfo,ninetyfiveErrors,signifMap,Deviance,residDevTestp] = ...
%                         FitLFPRetinoModel_Loglog(Data,xaxis,yaxis,numRepeats);
%                         MakePlots(finalParameters,meanResponse,xaxis,yaxis,stimLen,Radius,centerVals,numStimuli,numChans,jj,numFiles,h,ConditionNumber);
%                end
               
               dailyParameters{jj} = finalParameters;
               mapData{jj} = Data;
               parameterCI{jj} = ninetyfiveErrors;
               residualDevianceTest_pVal{jj} = residDevTestp;
               fisherInformation{jj} = fisherInfo;
               fullDeviance{jj} = Deviance;
               
               srpResponse = zeros(numChans,srp_reps,stimLen);
               meanVEP = zeros(numChans,stimLen);
               signifVEP = zeros(numChans,1);
               for kk=1:numChans
                   stimStrobes = strobeTimes(svStrobed == srp_word);
                   for ll=1:srp_reps
                       stimOnset = stimStrobes(ll);
                       [~,index] = min(abs(timeStamps-stimOnset));
                       temp = ChanData(index:index+stimLen-1,currentChannel(kk));
                       srpResponse(kk,ll,:) = temp;
                   end
                   allVEPs = squeeze(srpResponse(kk,:,:));
                   meanVEP(kk,:) = mean(allVEPs,1);
                   [~,minLatency] = min(meanVEP(kk,:));
                   srpTrialVEPs{jj,kk} = allVEPs;
                   minWin = (minLatency-30):(minLatency+20);
                   maxWin = (minLatency+30):(minLatency+130);

                   subplot(numFiles,numChans*2,kk+numChans+(jj-1)*numChans*2);plot(meanVEP(kk,:));axis([0 stimLen -400 200]);
                   xlabel('Time from Flip/Flop (ms)');ylabel('LFP Mag (\muVolts)');
                   title(sprintf('SRP VEP: Day %d',jj));
               end
               srpSize{jj} = max(srpResponse(:,:,maxWin),[],3)-min(srpResponse(:,:,minWin),[],3);
               srpVEP{jj} = meanVEP;
               significantVEP{jj} = signifVEP;
           elseif ConditionNumber == 3
               srpResponse = zeros(numChans,srp_reps,stimLen);
               meanVEP = zeros(numChans,stimLen);
               signifVEP = zeros(numChans,1);

               for kk=1:numChans
                   stimStrobes = strobeTimes(svStrobed == srp_word);
                   for ll=1:srp_reps
                       stimOnset = stimStrobes(ll);
                       [~,index] = min(abs(timeStamps-stimOnset));
                       temp = ChanData(index:index+stimLen-1,currentChannel(kk));
                       srpResponse(kk,ll,:) = temp;
                   end
                   allVEPs = squeeze(srpResponse(kk,:,:));
                   meanVEP(kk,:) = mean(allVEPs,1);
                   [~,minLatency] = min(meanVEP(kk,:));
                   srpTrialVEPs{jj,kk} = allVEPs;
                   minWin = (minLatency-30):(minLatency+20);
                   maxWin = (minLatency+30):(minLatency+130);
  
                   
                   subplot(numFiles,numChans,kk+(jj-1)*numChans);plot(meanVEP(kk,:));axis([0 stimLen -400 200]);
                   xlabel('Time from Flip/Flop (ms)');ylabel('LFP Mag (\muVolts)');
                   title(sprintf('SRP VEP: Day %d',jj));
               end
               srpSize{jj} = max(srpResponse(:,:,maxWin),[],3)-min(srpResponse(:,:,minWin),[],3);
               srpVEP{jj} = meanVEP;
               significantVEP{jj} = signifVEP;
           end
       end
%        savefig(h,sprintf('MappingEffectsResults_%d.fig',Animals(ii)));
        clear h;
   end
   filename = sprintf('MappingEffectsGammaTestResults_%d.mat',Animals(ii));
   if ConditionNumber == 1
       save(filename,'dailyParameters','parameterCI','fisherInformation',...
            'ConditionNumber','numFiles','w_pixels',...
            'h_pixels','mapData','residualDevianceTest_pVal','pix_to_degree',...
            'fullDeviance','currentChannel','mappingVEP');
   elseif ConditionNumber == 2
        save(filename,'dailyParameters','parameterCI','fisherInformation',...
            'srpSize','srpVEP','ConditionNumber','numFiles','w_pixels',...
            'h_pixels','mapData','residualDevianceTest_pVal','significantVEP',...
            'fullDeviance','pix_to_degree','currentChannel','mappingVEP','srpTrialVEPs');
   elseif ConditionNumber == 3
       save(filename,...
            'srpSize','srpVEP','ConditionNumber','numFiles','significantVEP','srpTrialVEPs');
   end
   fprintf('Done with animal %d\n\n',Animals(ii));pause(1);
end


end

function [ChanData,timeStamps,tsevs,svStrobed,numChans,sampleFreq] = ExtractSignal(EphysFileName)
    % Extract LFP signals from allad, filter, get timestamps
    % read in the .plx file

    if exist(strcat(EphysFileName(1:end-4),'.mat'),'file') ~= 2
        readall(EphysFileName);pause(0.1);
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
        n = 100;
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
    
end

function [] = MakePlots(finalParameters,meanResponse,xaxis,yaxis,stimLen,Radius,centerVals,numStimuli,numChans,Day,numFiles,h,ConditionNumber)
    figure(h);
    Radius = round(Radius);
    xPos = unique(centerVals(:,1));
    yPos = unique(centerVals(:,2));
    
    xDiff = mean(diff(xPos));
    yDiff = mean(diff(yPos));
    xconv = stimLen/xDiff; 
    yconv = 800/yDiff; % for height of the stimulus
    
    if ConditionNumber == 1
        numRows = numChans;
    elseif ConditionNumber == 2
        numRows = numChans*2;
    end
    
    for ii=1:numChans
        subplot(numFiles,numRows,ii+(Day-1)*numRows);axis([0 max(xaxis) 0 max(yaxis)]);
        if ii==1
            title(sprintf('VEP Retinotopy, Day %d',Day));
            xlabel('Horizontal Screen Position (pixels)');ylabel('Vertical Screen Position');
        end
        hold on;
%         figure();axis([0 max(xaxis) 0 max(yaxis)]);
        for jj=1:numStimuli
            tempx = centerVals(jj,1);
            tempy = centerVals(jj,2);
            plot(((1:1:stimLen)./xconv+centerVals(jj,1)-0.5*xDiff),...
                (squeeze(meanResponse(ii,jj,:))'./yconv+centerVals(jj,2)),'k','LineWidth',2);
            hold on;
            plot((tempx-Radius)*ones(Radius*2+1,1),(tempy-Radius):(tempy+Radius),'k','LineWidth',2);
            plot((tempx+Radius)*ones(Radius*2+1,1),(tempy-Radius):(tempy+Radius),'k','LineWidth',2);
            plot((tempx-Radius):(tempx+Radius),(tempy-Radius)*ones(Radius*2+1,1),'k','LineWidth',2);
            plot((tempx-Radius):(tempx+Radius),(tempy+Radius)*ones(Radius*2+1,1),'k','LineWidth',2);
            hold on;

        end
        
        finalIm = zeros(length(xaxis),length(yaxis));
        parameterVec = finalParameters(ii,:);
        b = [parameterVec(1),parameterVec(4),parameterVec(5),parameterVec(6)];
        for jj=1:length(xaxis)
            for kk=1:length(yaxis)
                distX = xaxis(jj)-parameterVec(2);
                distY = yaxis(kk)-parameterVec(3);
                
                finalIm(jj,kk) = b(1)*exp(-(distX.^2)./(2*b(2)*b(2))-...
                    (distY.^2)./(2*b(3)*b(3)))+b(4);
            end
        end
        imagesc(xaxis,yaxis,finalIm','AlphaData',0.5);set(gca,'YDir','normal');w=colorbar;
        caxis([b(4) b(4)+150]);
        ylabel(w,'Log Mean VEP Magnitude (\muV)');colormap('jet');
        hold off;
    end
end


