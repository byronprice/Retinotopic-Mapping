function [ ] = MappingEffects_GLM(Animals,Channels)
%MappingEffects_GLM.m
%  Analyze data from mapping effects experiment, see MappingEffects.m
%  Three experimental conditions:
%   1) retinotopic mapping, followed by fake SRP (grey screen)
%   2) retinotopic mapping, followed by SRP
%   3) fake mapping (grey screen), followed by SRP
%
% Created: 2017/04/14
%  Byron Price
% Updated: 2017/04/14
% By: Byron Price

cd('~/CloudStation/ByronExp/MappingEffects');

for ii=1:length(Animals)
   dataFiles = dir(sprintf('MappingEffectsData*%d.plx',Animals(ii)));
   stimFiles = dir(sprintf('MappingEffectsStim*%d.mat',Animals(ii)));
   
   numFiles = length(dataFiles);
   fprintf('Running animal %d\n\n',Animals(ii));
   
   linkFun = -2;
   
   pix_to_degree = cell(numFiles,1);
   currentChannel = Channels{ii};
   numChans = length(currentChannel);
   B = cell(numFiles,numChans);
   DEV = cell(numFiles,numChans,2);
   SE = cell(numFiles,numChans);
   DEV_RESIDUALS = cell(numFiles,numChans);
   F_TEST = cell(numFiles,numChans);
   VEP_MAGNITUDES = cell(numFiles,numChans);
   STIMULUS_LOCATIONS = cell(numFiles,numChans);
   ALL_VEPS = cell(numFiles,1);
   
   % could add back check for "goodChannels"
   if numChans ~= 0
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
               meanVEPMagnitude = zeros(numChans,numStimuli);
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
                       [~,minLatency] = min(squeeze(meanResponse(kk,ll,:)));
                       minLatency = max(min(minLatency,200),40);
                       minWin = minLatency-30:minLatency+30;
                       maxWin = minLatency+30:minLatency+200;
                       meanVEPMagnitude(kk,ll) = abs(max(squeeze(meanResponse(kk,ll,maxWin)))-min(squeeze(meanResponse(kk,ll,minWin))));
                   end
                   ALL_VEPS{jj} = squeeze(Response);
                   X = [centerVals(:,1),centerVals(:,2),centerVals(:,1).^2,centerVals(:,2).^2];
                   [b,dev,stats] = glmfit(X,meanVEPMagnitude(kk,:)','normal','link',linkFun);
                   B{jj,kk} = b;
                   DEV{jj,kk,1} = dev;
                   SE{jj,kk} = stats.se;
                   DEV_RESIDUALS{jj,kk} = stats.residd;
                   
                   [~,nullDev,~] = glmfit(ones(numStimuli,1),meanVEPMagnitude(kk,:)','normal','link',linkFun,'constant','off');
                   dev = real(dev);nullDev = real(nullDev);b = real(b);
                   fStat = ((nullDev-dev)/(length(b)-1))/(dev/(numStimuli-length(b)));
                   F_TEST{jj,kk} = fcdf(fStat,length(b)-1,numStimuli-length(b),'upper');
                   
                   VEP_MAGNITUDES{jj,kk} = meanVEPMagnitude(kk,:)';
                   STIMULUS_LOCATIONS{jj,kk} = centerVals;
                   
                   DEV{jj,kk,2} = nullDev;
                   
                   Radius = round(Radius);
                   xPos = unique(centerVals(:,1));
                   yPos = unique(centerVals(:,2));
                   
                   xDiff = mean(diff(xPos));
                   yDiff = mean(diff(yPos));
                   xconv = stimLen/xDiff;
                   yconv = 800/yDiff; % for height of the stimulus
                   figure();
                   for ww=1:numStimuli
                       tempx = centerVals(ww,1);
                       tempy = centerVals(ww,2);
                       plot(((1:1:stimLen)./xconv+centerVals(ww,1)-0.5*xDiff),...
                           (squeeze(meanResponse(kk,ww,:))'./yconv+centerVals(ww,2)),'k','LineWidth',2);
                       hold on;
                       plot((tempx-Radius)*ones(Radius*2+1,1),(tempy-Radius):(tempy+Radius),'k','LineWidth',2);
                       plot((tempx+Radius)*ones(Radius*2+1,1),(tempy-Radius):(tempy+Radius),'k','LineWidth',2);
                       plot((tempx-Radius):(tempx+Radius),(tempy-Radius)*ones(Radius*2+1,1),'k','LineWidth',2);
                       plot((tempx-Radius):(tempx+Radius),(tempy+Radius)*ones(Radius*2+1,1),'k','LineWidth',2);
                       hold on;
                       
                   end
                   
                   xPos = 1:w_pixels;yPos = 1:h_pixels;
                   [XPOS,YPOS] = meshgrid(xPos,yPos);
                   tempX = [XPOS(:),YPOS(:),XPOS(:).^2,YPOS(:).^2];
                   nx = length(xPos);ny = length(yPos);
                   yhat = glmval(b,tempX,linkFun);
                   imagesc(xPos,yPos,reshape(yhat,[ny,nx]),'AlphaData',0.5);set(gca,'YDir','normal');w=colorbar;
                   ylabel(w,'Mean VEP Magnitude (\muV)');colormap('jet');
                   title(sprintf('F test p-value: %3.2e',F_TEST{jj,kk}));
                   axis([0 w_pixels 0 h_pixels]);
                   hold off;
               end

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
               meanVEPMagnitude = zeros(numChans,numStimuli);
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
                       [~,minLatency] = min(squeeze(meanResponse(kk,ll,:)));
                       minLatency = max(min(minLatency,200),40);
                       minWin = minLatency-30:minLatency+30;
                       maxWin = minLatency+30:minLatency+200;
                       meanVEPMagnitude(kk,ll) = abs(max(squeeze(meanResponse(kk,ll,maxWin)))-min(squeeze(meanResponse(kk,ll,minWin))));
                   end
                   ALL_VEPS{jj} = squeeze(Response);
                   X = [centerVals(:,1),centerVals(:,2),centerVals(:,1).^2,centerVals(:,2).^2];
                   [b,dev,stats] = glmfit(X,meanVEPMagnitude(kk,:)','normal','link',linkFun);
                   B{jj,kk} = b;
                   DEV{jj,kk,1} = dev;
                   SE{jj,kk} = stats.se;
                   DEV_RESIDUALS{jj,kk} = stats.residd;
                   
                   [~,nullDev,~] = glmfit(ones(numStimuli,1),meanVEPMagnitude(kk,:)','normal','link',linkFun,'constant','off');
                   dev = real(dev);nullDev = real(nullDev);b = real(b);
                   fStat = ((nullDev-dev)/(length(b)-1))/(dev/(numStimuli-length(b)));
                   F_TEST{jj,kk} = fcdf(fStat,length(b)-1,numStimuli-length(b),'upper');
                   
                   DEV{jj,kk,2} = nullDev;
                   
                   VEP_MAGNITUDES{jj,kk} = meanVEPMagnitude(kk,:)';
                   STIMULUS_LOCATIONS{jj,kk} = centerVals;
                   
                   Radius = round(Radius);
                   xPos = unique(centerVals(:,1));
                   yPos = unique(centerVals(:,2));
                   
                   xDiff = mean(diff(xPos));
                   yDiff = mean(diff(yPos));
                   xconv = stimLen/xDiff;
                   yconv = 800/yDiff; % for height of the stimulus
                   
                   figure();
                   for ww=1:numStimuli
                       tempx = centerVals(ww,1);
                       tempy = centerVals(ww,2);
                       plot(((1:1:stimLen)./xconv+centerVals(ww,1)-0.5*xDiff),...
                           (squeeze(meanResponse(kk,ww,:))'./yconv+centerVals(ww,2)),'k','LineWidth',2);
                       hold on;
                       plot((tempx-Radius)*ones(Radius*2+1,1),(tempy-Radius):(tempy+Radius),'k','LineWidth',2);
                       plot((tempx+Radius)*ones(Radius*2+1,1),(tempy-Radius):(tempy+Radius),'k','LineWidth',2);
                       plot((tempx-Radius):(tempx+Radius),(tempy-Radius)*ones(Radius*2+1,1),'k','LineWidth',2);
                       plot((tempx-Radius):(tempx+Radius),(tempy+Radius)*ones(Radius*2+1,1),'k','LineWidth',2);
                       hold on;
                       
                   end
                   
                   xPos = 1:w_pixels;yPos = 1:h_pixels;
                   [XPOS,YPOS] = meshgrid(xPos,yPos);
                   tempX = [XPOS(:),YPOS(:),XPOS(:).^2,YPOS(:).^2];
                   nx = length(xPos);ny = length(yPos);
                   yhat = glmval(b,tempX,linkFun);
                   imagesc(xPos,yPos,reshape(yhat,[ny,nx]),'AlphaData',0.5);set(gca,'YDir','normal');w=colorbar;
                   ylabel(w,'Mean VEP Magnitude (\muV)');colormap('jet');
                   title(sprintf('F test p-value: %3.2e',F_TEST{jj,kk}));
                   axis([0 w_pixels 0 h_pixels]);
                   hold off;
               end
              
           
           end
       end
%        savefig(h,sprintf('MappingEffectsResults_%d.fig',Animals(ii)));
   end
   
   filename = sprintf('MappingEffectsGLMResults_%d.mat',Animals(ii));
   if ConditionNumber == 1 || ConditionNumber == 2
       save(filename,'ConditionNumber','pix_to_degree','w_pixels','h_pixels','currentChannel','numChans',...
           'B','DEV','SE','DEV_RESIDUALS','F_TEST','VEP_MAGNITUDES','STIMULUS_LOCATIONS','ALL_VEPS');
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
