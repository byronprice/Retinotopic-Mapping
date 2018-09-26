function [] = MapRetinotopyCallawayOE(AnimalName,Date)
% MapRetinotopyCallawayOE.m
%
%  Will take data from a retinotopic mapping experiment and extract the
%   retinotopy of the LFP recording electrode. The stimulus used here is
%   the periodic drifting bar with flashing counter-phase checkerboard from
%   the file RetinotopyCallaway.m
%     For use with open ephys recording system
%
%INPUT: AnimalName - unique identifier for the animal as a number, e.g.
%            12345
%       Date - date of the experiment input as a number yearMonthDay, 
%            e.g. 20160525
%OUTPUT: plots
%
% Created: 2016/05/31, 24 Cummington, Boston
%  Byron Price
% Updated: 2018/09/25
%  By: Byron Price

% read in the .mat compiled data file, from open ephys
StimulusFileName = strcat('RetinoCallStim',num2str(Date),'_',num2str(AnimalName),'.mat');
eyeTrackFileName = strcat('RetinoCall_',num2str(Date),'-',num2str(AnimalName),'-MLP.mat');

newDate = num2str(Date);newDate = [newDate(1:4),'-',newDate(5:6),'-',newDate(7:8)];
EphysFileName = strcat('CompiledData_RetinoCall_',newDate,'*',num2str(AnimalName),'*');
file = dir(EphysFileName);

EphysFileName = file(1).name;
load(EphysFileName,'auxData','events','eventTimes','lowpassData','lowpassTimes',...
    'lpFs','numChans');numChans = 1;

load(StimulusFileName);
load(eyeTrackFileName,'blink','flagged','Fs','pupilRotation');

badEyeTrack = logical(blink+flagged);
medianEyePos = median(pupilRotation(~badEyeTrack,:),1);

pupilRotation(:,1) = pupilRotation(:,1)-medianEyePos(1); %left-right eye movements
pupilRotation(:,2) = -(pupilRotation(:,2)-medianEyePos(2));% up-down eye movements

pupilRotation(badEyeTrack,:) = NaN;

driftSpeed = stimParams.driftSpeed; % units of degrees per second
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
ifi = stimParams.ifi;
mmPerPixel = stimParams.mmPerPixel;
DistToScreen = stimParams.DistToScreen;
theta = stimParams.theta;
phi = stimParams.phi;

stimulationFrequency = 1/checkRefresh;

sampleFreq = lpFs;

stimLen = round(driftTime*sampleFreq);
stimLenDS = round(driftTime*Fs);

downsampleFactor = round(mean(stimLen./stimLenDS));

stimLenDS = ceil(stimLen./downsampleFactor);

vidInds = find(auxData(:,2));

% COLLECT DATA IN THE PRESENCE OF VISUAL STIMULI
Response = cell(numChans,4);
for ii=1:numChans
    Response{ii,1} = zeros(2*reps,stimLen(1));
    Response{ii,2} = zeros(2*reps,stimLen(2));
    Response{ii,3} = zeros(2*reps,stimLenDS(1));
    Response{ii,4} = zeros(2*reps,stimLenDS(2));
    
    horzCount = 1;
    vertCount = 1;
    for jj=1:numDirs
        strobeNum = jj;
        currentStrobeTimes = eventTimes(events==strobeNum);
        for kk=1:reps
           [~,onsetTime] = min(abs(currentStrobeTimes(kk)+ifi-lowpassTimes));
           [~,vidOn] = min(abs(onsetTime-vidInds));
           
           if strcmp(DirNames{jj},'Left') == 1 || strcmp(DirNames{jj},'Right') == 1
               offsetTime = onsetTime+stimLen(1)-1;
               vidOff = vidOn+stimLenDS(1)-1;
           elseif strcmp(DirNames{jj},'Down') == 1 || strcmp(DirNames{jj},'Up') == 1
               offsetTime = onsetTime+stimLen(2)-1;
               vidOff = vidOn+stimLenDS(2)-1;
           end
           
           if strcmp(DirNames{jj},'Left') == 1 || strcmp(DirNames{jj},'Down') == 1
              tempLFP = flipud(lowpassData(onsetTime:offsetTime,ii));
              tempEye = flipud(pupilRotation(vidOn:vidOff,:));
           else
              tempLFP = lowpassData(onsetTime:offsetTime,ii);
              tempEye = pupilRotation(vidOn:vidOff,:);
           end
           
           if strcmp(DirNames{jj},'Left') == 1 || strcmp(DirNames{jj},'Right') == 1
                Response{ii,1}(horzCount,:) = tempLFP;
                Response{ii,3}(horzCount,:) = tempEye(:,1);
                horzCount = horzCount+1;
           elseif strcmp(DirNames{jj},'Down') == 1 || strcmp(DirNames{jj},'Up') == 1
                Response{ii,2}(vertCount,:) = tempLFP;
                Response{ii,4}(vertCount,:) = tempEye(:,2);
                vertCount = vertCount+1;
           end
        end
    end
end

%time = linspace(0,driftTime,stimLen);
horzPosition = linspace(min(phi),max(phi),stimLen(1)); % azimuth
vertPosition = linspace(min(theta),max(theta),stimLen(2)); % altitude

waveletSize = 50; % 
kernelLen = round(waveletSize*checkRefresh*sampleFreq);
if mod(kernelLen,2) == 0
    kernelLen = kernelLen+1;
end
x = linspace(-waveletSize/2*checkRefresh,waveletSize/2*checkRefresh,kernelLen);

stdGauss = (waveletSize/2)*checkRefresh/4;
gaussKernel = exp(-(x.*x)./(2*stdGauss*stdGauss));
kernel = exp(-2*pi*x*1i*stimulationFrequency).*gaussKernel;

noiseFreqs = [stimulationFrequency-5,stimulationFrequency-2.5,stimulationFrequency-2,stimulationFrequency-1.5,...
    stimulationFrequency+1.5,stimulationFrequency+2,stimulationFrequency+2.5,stimulationFrequency+5];
noiseKernels = zeros(length(noiseFreqs),length(kernel));
for ii=1:length(noiseFreqs)
    noiseKernels(ii,:) = exp(-2*pi*x*1i*noiseFreqs(ii)).*gaussKernel;
end

transformBaseline = zeros(numChans,1);
for ii=1:numChans
    transformBaseline(ii) = 1;
end

DirNames = cell(2,1);DirNames{1} = 'Horizontal Sweep';DirNames{2} = 'Vertical Sweep';
transformResponse = cell(numChans,2);
Results = struct('b',{cell(numChans,2)},'se',{cell(numChans,2)},...
    'ScreenPos',{cell(numChans,2)},'Center',{cell(numChans,2)},...
    'FWHM',{cell(numChans,2)},'Dev',{cell(numChans,2)},'Data',{cell(numChans,2)});

for ii=1:numChans
    for jj=1:2
        transformResponse{ii,jj} = zeros(2*reps,stimLenDS(jj));
        for kk=1:2*reps
            data = Response{ii,jj}(kk,:);
            convData = conv(data,kernel,'same');
            %             convData = convData(1:length(data));
            convData = sqrt(convData.*conj(convData));
            
            noiseConvData = zeros(size(convData));
            for ll=1:length(noiseFreqs)
                temp = conv(data,noiseKernels(ll,:),'same');
                temp = sqrt(temp.*conj(temp));
                noiseConvData = noiseConvData+temp./length(noiseFreqs);
            end
            
            transformResponse{ii,jj}(kk,:) = convData(1:downsampleFactor:end)./...
                noiseConvData(1:downsampleFactor:end)-transformBaseline(ii);
        end
    end
end

radianPerPixel = (0.25:0.05:10).*pi./180;

deviances = zeros(length(radianPerPixel),1);

for zz=1:length(radianPerPixel)
    for ii=1:numChans
        Y = [];
        POS = [];
        for jj=1:2
            position = zeros(2*reps,stimLenDS(jj));
            for kk=1:2*reps
                
                if jj==1
                    tmp = Response{ii,jj+2}(kk,:);
                    position(kk,:) = horzPosition(1:downsampleFactor:end)+tmp.*radianPerPixel(zz);
                elseif jj==2
                    tmp = Response{ii,jj+2}(kk,:);
                    position(kk,:) = vertPosition(1:downsampleFactor:end)+tmp.*radianPerPixel(zz);
                end
            end
            y = transformResponse{ii,jj};y=y(:);
            
            Y = [Y;y];
            if jj==1
                POS = [POS;position(:),zeros(length(position(:)),1)];
            elseif jj==2
                POS = [POS;zeros(length(position(:)),1),position(:)];
            end
            
        end
        inds = find(~isnan(sum(POS,2)));
        Y = Y(inds);POS = POS(inds,:);
        
        design = [ones(length(Y),1),POS(:,1),POS(:,1).^2,POS(:,2),POS(:,2).^2];
        [~,mainDev,~] = glmfit(design,Y,'normal','link','log','constant','off');
        
%         fprintf('Deviance: %3.2f\n',real(mainDev));
        deviances(zz) = real(mainDev);
    end
end

[~,ind] = min(deviances);
radianPerPixel = radianPerPixel(ind);
fprintf('Degrees per pixel: %3.2f\n',radianPerPixel*180/pi);

for ii=1:numChans
    Y = [];
    POS = [];
    for jj=1:2
        position = zeros(2*reps,stimLenDS(jj));
        for kk=1:2*reps
            
            if jj==1
                tmp = Response{ii,jj+2}(kk,:);
                position(kk,:) = horzPosition(1:downsampleFactor:end)+tmp.*radianPerPixel;
            elseif jj==2
                tmp = Response{ii,jj+2}(kk,:);
                position(kk,:) = vertPosition(1:downsampleFactor:end)+tmp.*radianPerPixel;
            end
        end
        y = transformResponse{ii,jj};y=y(:);
        
        Y = [Y;y];
        if jj==1
            POS = [POS;position(:),zeros(length(position(:)),1)];
        elseif jj==2
            POS = [POS;zeros(length(position(:)),1),position(:)]; 
        end
        
    end
    inds = find(~isnan(sum(POS,2)));
    Y = Y(inds);POS = POS(inds,:);
    Results.Data{ii,1} = Y;
    Results.Data{ii,2} = POS;
    design = [ones(length(Y),1),POS(:,1),POS(:,1).^2,POS(:,2),POS(:,2).^2];
    [b,mainDev,stats] = glmfit(design,Y,'normal','link','log','constant','off');
    
    b = real(b);
    standardError = real(stats.se);
    
    Results.b{ii,1} = b;
    Results.se{ii,1} = standardError;
    
    restrictDesign = ones(length(Y),1);
    [~,restrictDev,~] = glmfit(restrictDesign,Y,'normal','constant','off');
    
    Results.Dev{ii,1} = real([mainDev,restrictDev]);
    
    Results.ScreenPos{ii,1} = unique(POS(:,1));
    Results.Center{ii,1} = (-b(2)/(2*b(3)));
    Results.FWHM{ii,1} = 2*sqrt(-log(2)/b(3));
    
    Results.ScreenPos{ii,2} = unique(POS(:,2));
    Results.Center{ii,2} = (-b(4)/(2*b(5)));
    Results.FWHM{ii,2} = 2*sqrt(-log(2)/b(5));
    
%     for jj=1:2
%         forDisplayDesign = [ones(size(position,2),1),position(1,:)',position(1,:)'.*position(1,:)'];
%         
%         yhat = exp(forDisplayDesign*b);
%         
%         temp = position';y = reshape(y,[2*reps,dsStimLen(jj)])';y = y(:);
%         temp = temp(:);
%         midPoint = round(length(y)/2);
%         figure();plot(temp(1:midPoint),y(1:midPoint),'.b');
%         hold on;
%         plot(temp(midPoint+1:end),y(midPoint+1:end),'.b');hold on;
%         plot(position(1,:)',yhat,'c','LineWidth',5)
%         title(sprintf('Chan: %d - %s',ii,DirNames{jj}));
%     end
end

fprintf('X Center: %3.2f\n',Results.Center{ii,1}.*180/pi);
fprintf('Y Center: %3.2f\n',Results.Center{ii,2}.*180/pi);

figure();
for ii=1:numChans
   xPos = Results.ScreenPos{ii,1};
   yPos = Results.ScreenPos{ii,2};
   b = Results.b{ii,1};

   POS = [xPos,zeros(length(xPos),1);zeros(length(yPos),1),yPos];
   
   Design = [ones(length(POS),1),POS(:,1),POS(:,1).^2,POS(:,2),POS(:,2).^2];
   
   muHorz = exp(Design(1:length(xPos),:)*b);
   muVert = exp(Design(length(xPos)+1:end,:)*b);
   
   muHorz = repmat(muHorz',[length(muVert),1]);
   muVert = repmat(muVert,[1,length(muHorz)]);
   
   finalIm = muHorz.*muVert;
   subplot(numChans,1,ii);
   imagesc(linspace(xPos(1),xPos(end),length(muHorz)).*180/pi,linspace(yPos(1),yPos(end),length(muVert)).*180/pi,finalIm);
   set(gca,'YDir','normal');colormap('jet');
   title(sprintf('LFP Retinotopy: Chan %d, Animal %d',ii,AnimalName));
   xlabel('Azimuth (degrees)');
   ylabel('Altitude (degrees)');
end

fileName = sprintf('RetinoCallResults%d_%d.mat',Date,AnimalName);
save(fileName,'transformResponse','Results','DirNames','w_pixels','h_pixels',...
    'stimulationFrequency','waveletSize','kernel','Response','stimLen',...
    'stimLenDS','downsampleFactor','centerPos','mmPerPixel','ifi','DistToScreen',...
    'radianPerPixel');

end
