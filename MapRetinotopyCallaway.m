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
%OUTPUT: plots
%
% Created: 2016/05/31, 24 Cummington, Boston
%  Byron Price
% Updated: 2017/10/22
%  By: Byron Price

% read in the .plx file
EphysFileName = strcat('RetinoCallData',num2str(Date),'_',num2str(AnimalName));

if exist(strcat(EphysFileName,'.mat'),'file') ~= 2
    readall(strcat(EphysFileName,'.plx'));
end

StimulusFileName = strcat('RetinoCallStim',num2str(Date),'_',num2str(AnimalName),'.mat');
EphysFileName = strcat(EphysFileName,'.mat');
load(EphysFileName,'allad','adfreq','tsevs','svStrobed','SlowPeakV',...
    'SlowADResBits','adgains','adfreqs','allts')
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
    ChanData(:,ii) = voltage;
%     n = 20;
%     lowpass = 100/(sampleFreq/2); % fraction of Nyquist frequency
%     blo = fir1(n,lowpass,'low',hamming(n+1));
%     temp = filter(blo,1,voltage);
%     
%     notch = 60/(sampleFreq/2);
%     bw = notch/n;
%     [b,a] = iirnotch(notch,bw);
%     ChanData(:,ii) = filter(b,a,temp);
end


timeStamps = 0:1/sampleFreq:dataLength/sampleFreq-1/sampleFreq;

if length(timeStamps) ~= dataLength
    display('Error: Review allad cell array and timing')
    return;
end
strobeTimes = tsevs{1,strobeStart};
stimLen = round(driftTime*sampleFreq);

% COLLECT DATA IN THE PRESENCE OF VISUAL STIMULI
Response = cell(numChans,3);

for ii=1:numChans
    Response{ii,1} = zeros(2*reps,stimLen(1));
    Response{ii,2} = zeros(2*reps,stimLen(2));
    
    greyStrobes = strobeTimes(svStrobed==0);
    numGrey = length(greyStrobes);
    Response{ii,3} = zeros(2*numGrey,stimLen(1));
    
    count = 1;
    for jj=1:numGrey
       onsetTime = round(greyStrobes(jj)*sampleFreq); 
       offsetTime = onsetTime+stimLen(1)-1;
       Response{ii,3}(count,:) = ChanData(onsetTime:offsetTime,ii);
       count = count+1;
       onsetTime = offsetTime+1;
       offsetTime = onsetTime+stimLen(1)-1;
       Response{ii,3}(count,:) = ChanData(onsetTime:offsetTime,ii);
       count = count+1;
    end
    
    horzCount = 1;
    vertCount = 1;
    for jj=1:numDirs
        strobeNum = jj;
        currentStrobeTimes = strobeTimes(svStrobed==strobeNum);
        for kk=1:reps
           onsetTime = round(currentStrobeTimes(kk)*sampleFreq);
           
           if strcmp(DirNames{jj},'Left') == 1 || strcmp(DirNames{jj},'Right') == 1
               offsetTime = onsetTime+stimLen(1)-1;
           elseif strcmp(DirNames{jj},'Down') == 1 || strcmp(DirNames{jj},'Up') == 1
               offsetTime = onsetTime+stimLen(2)-1;
           end
           
           if strcmp(DirNames{jj},'Left') == 1 || strcmp(DirNames{jj},'Down') == 1
              tempLFP = flipud(ChanData(onsetTime:offsetTime,ii));
           else
              tempLFP = ChanData(onsetTime:offsetTime,ii);
           end
           
           if strcmp(DirNames{jj},'Left') == 1 || strcmp(DirNames{jj},'Right') == 1
                Response{ii,1}(horzCount,:) = tempLFP;
                horzCount = horzCount+1;
           elseif strcmp(DirNames{jj},'Down') == 1 || strcmp(DirNames{jj},'Up') == 1
                Response{ii,2}(vertCount,:) = tempLFP;
                vertCount = vertCount+1;
           end
        end
    end
end

%time = linspace(0,driftTime,stimLen);
horzPosition = linspace(0,1,stimLen(1)); %w_pixels
vertPosition = linspace(0,1,stimLen(2)); %h_pixels

waveletSize = 20;
x = linspace(-waveletSize/2*checkRefresh,waveletSize/2*checkRefresh,round(waveletSize*checkRefresh*sampleFreq));
kernel = exp(-2*pi*x*1i*stimulationFrequency);
transformBaseline = zeros(numChans,1);
for ii=1:numChans
    numGrey = size(Response{ii,3},1);
    temp = zeros(numGrey,stimLen(1));
    for jj=1:numGrey
        data = Response{ii,3}(jj,:);
        data = conv(data,kernel,'same');
        data = sqrt(data.*conj(data));
        temp(jj,:) = data;
    end
    transformBaseline(ii) = mean(temp(:));
end

DirNames = cell(2,1);DirNames{1} = 'Horizontal Sweep';DirNames{2} = 'Vertical Sweep';
transformResponse = cell(numChans,2);
Results = struct('b',{cell(numChans,2)},'se',{cell(numChans,2)},...
    'F',{cell(numChans,2)},'ScreenPos',{cell(numChans,2)},'Center',{cell(numChans,2)},...
    'FWHM',{cell(numChans,2)});

downsampleFactor = 10;
dsStimLen = ceil(stimLen/downsampleFactor);
for ii=1:numChans
    for jj=1:2
        transformResponse{ii,jj} = zeros(2*reps,dsStimLen(jj));
        position = zeros(2*reps,dsStimLen(jj));
        for kk=1:2*reps
            data = Response{ii,jj}(kk,:);
            data = conv(data,kernel,'same');
            data = sqrt(data.*conj(data));
            %figure();parcorr(data(1:downsampleFactor:end))
            transformResponse{ii,jj}(kk,:) = data(1:downsampleFactor:end)-transformBaseline(ii);
            
            if jj==1
                position(kk,:) = horzPosition(1:downsampleFactor:end);
            elseif jj==2
                position(kk,:) = vertPosition(1:downsampleFactor:end);
            end
        end
        y = transformResponse{ii,jj};y=y(:);

        effectiveN = length(y);
        
        design = [ones(length(y),1),position(:),position(:).*position(:)];
        [b,~,stats] = glmfit(design,y,'normal','link','log','constant','off');

        %[b,FI] = GetMLest(design,y,b);
        
        %asymptotVar = pinv((1/effectiveN).*FI);
        %standardError = sqrt(diag(asymptotVar));
        
        b = real(b);
        standardError = real(stats.se);
        
        Results.b{ii,jj} = b;
        Results.se{ii,jj} = standardError;
        
        mainDev = GetDeviance(design,y,b);
        
        % F-test
        restrictDesign = ones(length(y),1);
        [~,restrictDev,~] = glmfit(restrictDesign,y,'normal','constant','off');
        
        F = ((restrictDev-mainDev)/(length(b)-1))/(mainDev/(effectiveN-length(b)-1));
        Ftest_p = fcdf(F,length(b)-1,effectiveN-length(b),'upper');
        
        Results.F{ii,jj} = [F,Ftest_p,length(b)-1,effectiveN-length(b)];
        
        if jj==1
            Results.ScreenPos{ii,jj} = position(1,:)'.*w_pixels;
            Results.Center{ii,jj} = -b(2)/(2*b(3))*w_pixels;
            Results.FWHM{ii,jj} = 2*sqrt(-log(2)/b(3))*w_pixels;
        elseif jj==2
            Results.ScreenPos{ii,jj} = position(1,:)'.*h_pixels; 
            Results.Center{ii,jj} = -b(2)/(2*b(3))*h_pixels;
            Results.FWHM{ii,jj} = 2*sqrt(-log(2)/b(3))*h_pixels;
        end
        
        forDisplayDesign = [ones(size(position,2),1),position(1,:)',position(1,:)'.*position(1,:)'];

        yhat = exp(forDisplayDesign*b);
        
        temp = position';y = reshape(y,[2*reps,dsStimLen(jj)])';y = y(:);
        temp = temp(:);
        midPoint = round(length(y)/2);
        figure();plot(temp(1:midPoint),y(1:midPoint),'.b');
        hold on;
        plot(temp(midPoint+1:end),y(midPoint+1:end),'.b');hold on;
        plot(position(1,:)',yhat,'c','LineWidth',5)
        title(sprintf('Chan: %d - %s',ii,DirNames{jj}));
    end
end

for ii=1:numChans
   xPos = Results.ScreenPos{ii,1};
   yPos = Results.ScreenPos{ii,2};
   bHorz = Results.b{ii,1};
   bVert = Results.b{ii,2};
   
   horzDesign = [ones(length(xPos),1),xPos./w_pixels,(xPos./w_pixels).^2];
   vertDesign = [ones(length(xPos),1),yPos./h_pixels,(yPos./h_pixels).^2];
   
   muHorz = exp(horzDesign*bHorz);
   muVert =exp(vertDesign*bVert);
   
   muHorz = repmat(muHorz',[dsStimLen(1),1]);
   muVert = repmat(muVert,[1,dsStimLen(2)]);
   
   finalIm = muHorz.*muVert;
   figure();imagesc(linspace(0,w_pixels,dsStimLen(1)),linspace(0,h_pixels,dsStimLen(2)),finalIm);
   set(gca,'YDir','normal');colormap('jet');
   title(sprintf('LFP Retinotopy: Chan %d, Animal %d',ii,AnimalName));
   xlabel('Horizontal Screen Position (pixels)');
   ylabel('Vertical Screen Position (pixels)');
end

fileName = sprintf('RetinoCallResults%d_%d.mat',Date,AnimalName);
save(fileName,'transformResponse','Results','DirNames','w_pixels','h_pixels',...
    'stimulationFrequency','waveletSize','kernel','Response','stimLen');

% timebandwidth = 60; % approximate standard deviation in time is 
%                 % 0.5*sqrt(timebandwidth/2)
% for ii=1:numChans
%     for jj=1:2
%         figure();
%         for kk=1:2*reps
%             [wt,f] = cwt(Response{ii,jj}(kk,:),sampleFreq,'TimeBandwidth',timebandwidth);
%             wt = sqrt(wt.*conj(wt));
%             [~,ind] = min(abs(f-stimulationFrequency));
%             wt = wt(ind:ind,:);
%             wt = mean(wt,1);
%             baseline = quantile(wt,0.05);
%             subplot(2*reps,1,kk);plot(time,wt-baseline);
%         end
%     end
% end

% discrete or continuous wavelet transform, hilbert transform
%  hilbert(x), dwt, or cwt
% windowLen = floor(checkRefresh*sampleFreq);
% x = 0:windowLen;y = sin(2*pi*x/windowLen);
% z = Response{1,1}(1,:);
% w = conv(z,y);
% plot(w.*w);
end

function [b,FI] = GetMLest(design,y,b)
b = real(b);

maxIter = 3e3;
tolerance = 1e-6;

stepSize = 1e-2;
numParams = length(b);

currentParams = b;
currentLikelihood = GetNormalLikelihood(design,y,currentParams);

iter = 1;
difference = 1;
lineSteps = [stepSize,1e-6,1e-5,1e-4,1e-3,1e-1,0.5e-1,1,10];lineN = length(lineSteps);
lineLikelies = zeros(lineN,1);
while iter < maxIter && difference > tolerance
    difference = 0;
    for jj=1:numParams
       newB = currentParams;newB(jj) = currentParams(jj)+lineSteps(1);
       newLikelihood = GetNormalLikelihood(design,y,newB);
       
       gradient = (newLikelihood-currentLikelihood)./lineSteps(1);
       lineLikelies(1) = newLikelihood;
       for kk=2:lineN
          newB = currentParams;newB(jj) = currentParams(jj)+sign(gradient)*lineSteps(kk);
          lineLikelies(kk) = GetNormalLikelihood(design,y,newB);
       end
       
       [maxLikely,ind] = max(lineLikelies);
       if maxLikely > currentLikelihood
           difference = difference+maxLikely-currentLikelihood;
           currentLikelihood = maxLikely;
           currentParams(jj) = currentParams(jj)+sign(gradient)*lineSteps(ind);
       else
           difference = difference+1;
       end
       iter = iter+1;
    end
end
b = currentParams;
FI = GetFisherInfo(design,y,b);
end

function [loglikelihood] = GetNormalLikelihood(design,y,b)
n = length(y);

% log link in this case
sumsquares = sum((exp(design*b)-y).^2);
residualVar = (1/n)*sumsquares;
loglikelihood = -(n/2)*log(2*pi)-(n/2)*log(residualVar)-(1/(2*residualVar))*sumsquares;

end

function [deviance] = GetDeviance(design,y,b)

deviance = sum((exp(design*b)-y).^2);

end

function [FI] = GetFisherInfo(design,y,b)
stepSize = 1e-3;
numParams = length(b);

set = [1,1;1,-1;-1,1;-1,-1];
signs = [1,-1,-1,1];
FI = zeros(numParams,numParams);
for ii=1:numParams
    for jj=1:numParams
        temp = 0;
        for kk=1:4
           shifts = set(kk,:);
           newB = b;newB(ii) = b(ii)+shifts(1)*stepSize;
           newB(jj) = b(jj)+shifts(2)*stepSize;
           
           temp = temp+signs(kk)*GetNormalLikelihood(design,y,newB);
        end
        
        FI(ii,jj) = (1/(4*stepSize*stepSize))*(temp);
    end
end

end
