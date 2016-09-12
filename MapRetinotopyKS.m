function [] = MapRetinotopyKS(AnimalName,Date)
% MapRetinotopyKS.m
%
%  Will take data from a retinotopic mapping experiment and extract the
%   retinotopy of the LFP recording electrode using a KS test, which
%   compares the distribution of latency to minimum and latency to maximum
%   against a uniform distribution (under the null hypothesis that no VEP
%   occurred, the distribution of latency times should be uniform). 
%INPUT: AnimalName - unique identifier for the animal as a number, e.g.
%            12345
%       Date - date of the experiment, e.g. 20160525

%OUTPUT: saved files and figures with info regarding retinotopy of each
%         channel
%
% Created: 2016/08/29, 24 Cummington, Boston
%  Byron Price
% Updated: 2016/08/29
%  By: Byron Price

cd('~/CloudStation/ByronExp/Retino');

EphysFileName = sprintf('RetinoData%d_%d',Date,AnimalName); % no file identifier
    % because MyReadall does that for us
  
global centerVals Radius reps stimTime holdTime numStimuli w_pixels h_pixels ...
    DistToScreen numChans sampleFreq stimLen maxWin minWin; %#ok<*REDEF>

StimulusFileName = sprintf('RetinoStim%d_%d.mat',Date,AnimalName);
load(StimulusFileName)

centerVals = stimParams.centerVals;
Radius = stimParams.Radius;
reps = stimParams.reps;
stimTime = stimParams.stimTime;
holdTime = stimParams.holdTime;
numStimuli = stimParams.numStimuli;
w_pixels = stimParams.w_pixels;
h_pixels = stimParams.h_pixels;
DistToScreen = stimParams.DistToScreen;

reps = reps-1;
% convert from allad to ChanData by filtering
[ChanData,timeStamps,tsevs,svStrobed] = ExtractSignal(EphysFileName);

% get LFP response to each stimulus (the VEPs)
[Response,meanResponse,strobeTimes] = CollectVEPS(ChanData,timeStamps,tsevs,svStrobed);

% BOOTSTRAP FOR 95% CONFIDENCE INTERVALS OF STATISTIC IN ABSENCE OF VISUAL STIMULI
noStimLen = holdTime*sampleFreq-stimLen*2;

baseStats = struct;
baseStats.minci = zeros(numChans,2);
baseStats.min = zeros(numChans,1);
baseStats.minsem = zeros(numChans,1);
baseStats.maxci = zeros(numChans,2);
baseStats.max = zeros(numChans,1);
baseStats.maxsem = zeros(numChans,1);

pauseOnset = strobeTimes(svStrobed == 0);
nums = length(pauseOnset);
N = 2000;
alpha = 0.01;
for ii=1:numChans
    Tboot = zeros(N,2);
    for jj=1:N
        indeces = random('Discrete Uniform',noStimLen,[reps,1]);
        temp = zeros(reps,stimLen);
        num = random('Discrete Uniform',nums);
        [~,index] = min(abs(timeStamps-pauseOnset(num)));
        index = index+stimLen;
        for kk=1:reps
            temp(kk,:) = ChanData(index+indeces(kk):index+indeces(kk)+stimLen-1,ii);
        end
        [~,minlats] = min(temp(:,minWin),[],2);
        [~,maxlats] = max(temp(:,maxWin),[],2);
        Tboot(jj,1) = max(eCDF(minlats,minWin)-linspace(0,1,length(minWin))');
        Tboot(jj,2) = max(eCDF(maxlats,maxWin)-linspace(0,1,length(maxWin))');
    end
    baseStats.minci(ii,:) = [quantile(Tboot(:,1),alpha/2),quantile(Tboot(:,1),1-alpha/2)];
    baseStats.min(ii) = mean(Tboot(:,1));
    baseStats.minsem(ii) = std(Tboot(:,1));
    baseStats.maxci(ii,:) = [quantile(Tboot(:,2),alpha/2),quantile(Tboot(:,1),1-alpha/2)];
    baseStats.max(ii) = mean(Tboot(:,2));
    baseStats.maxsem(ii) = std(Tboot(:,2));
end

% KS TEST to determine which stimuli are significant
[significantStimuli,minLatency,maxLatency] = KStest(Response,meanResponse,baseStats);

% xPos = unique(centerVals(:,1));
% yPos = unique(centerVals(:,2));
% numX = length(xPos);
% numY = length(yPos);
% position = zeros(numStimuli,1);
% 
% mapping = zeros(numX,numY);
% count = 1;
% for ii=numY:-1:1
%     for jj=1:numX
%         mapping(jj,ii) = count;
%         count = count+1;
%     end
% end
% for ii=1:numStimuli
%     xInd = find(logical(centerVals(ii,1) == xPos));
%     yInd = find(logical(centerVals(ii,2) == yPos));
%     position(ii) = mapping(xInd,yInd);
% end
% 
% for ii=1:numChans
%     w(1+2*(ii-1)) = figure();
%     w(2+2*(ii-1)) = figure();
%     for jj=1:numStimuli
% %         AIC = zeros(1,4);
% %         GMModels = cell(1,4);
% %         options = statset('MaxIter',500);
% %         for kk=1:3
% %             GMModels{kk} = fitgmdist(squeeze(minLatency(ii,jj,:)),kk,'Options',options);
% %             AIC(kk)= GMModels{kk}.AIC;
% %         end
% %         [~,numComponents] = min(AIC);
%         figure(w(1+2*(ii-1)));
%         subplot(numY,numX,position(jj));histogram(squeeze(minLatency(ii,jj,:)));
%         legend(num2str(std(squeeze(minLatency(ii,jj,:)))));
%         figure(w(2+2*(ii-1)));
%         subplot(numY,numX,position(jj));histogram(squeeze(maxLatency(ii,jj,:)));
%     end
% end

% Bootstrap 95% confidence interval on the latency
alpha = 0.05;
N = 2000;

minInterval = zeros(numChans,numStimuli,2);
maxInterval = zeros(numChans,numStimuli,2);
for ii=1:numChans
    for jj=1:numStimuli
        Data = squeeze(Response(ii,jj,:,:));
        [~,~,minci] = Bootmin(Data,alpha,N,reps);
        [~,~,maxci] = Bootmax(Data,alpha,N,reps);
        minInterval(ii,jj,:) = minci;
        maxInterval(ii,jj,:) = maxci;
    end
end
% Calculate the center of mass of the receptive field
[centerMass] = GetReceptiveField(significantStimuli,AnimalName);

[stimVals,h] = MakePlots(significantStimuli,meanResponse,centerMass,AnimalName,minInterval,maxInterval); 

Channel = input('Type the channel that looks best (as a number, e.g. 1): ');
savefig(h,sprintf('RetinoMapKS%d_%d.fig',Date,AnimalName));
%print(h,'-depsc','filename');


MapParams = RetinoMapObj;
MapParams.numChans = numChans;
MapParams.centerVals = centerVals;
MapParams.significantStimuli = significantStimuli;
MapParams.centerMass = centerMass;
MapParams.stimVals = stimVals;
MapParams.Response = Response;
MapParams.meanResponse = meanResponse;
MapParams.Channel = Channel;
MapParams.minLatency = minLatency;
MapParams.maxLatency = maxLatency;

save(sprintf('RetinoMapKS%d_%d.mat',Date,AnimalName),'MapParams');

yesNo = input('Save as principal map parameter file? (y/n): ','s');

if strcmp(yesNo,'y') == 1
    save(sprintf('RetinoMapKS%d.mat',AnimalName),'MapParams');
end

% obj = gmdistribution(centerMass(Channel,1:2),squeeze(Sigma(Channel,:,:)));
% figure();
% h = ezcontour(@(x,y) pdf(obj,[x y]),[0 w_pixels,0 h_pixels]);
end

function [ChanData,timeStamps,tsevs,svStrobed] = ExtractSignal(EphysFileName)
    % Extract LFP signals from allad, filter, get timestamps
    global numChans sampleFreq;
    % read in the .plx file

    if exist(strcat(EphysFileName,'.mat'),'file') ~= 2
        MyReadall(EphysFileName);
    end

    EphysFileName = strcat(EphysFileName,'.mat');
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
        display('Error: Review allad cell array and timing')
        return;
    end
end

function [Response,meanResponse,strobeTimes] = CollectVEPS(ChanData,timeStamps,tsevs,svStrobed)
    global numChans numStimuli reps stimLen maxWin sampleFreq minWin;
    strobeStart = 33;
    strobeTimes = tsevs{1,strobeStart};
    stimLen = round(0.25*sampleFreq); % about 250 milliseconds
    maxWin = round(0.1*sampleFreq):stimLen;
    minWin = round(0.05*sampleFreq):round(0.2*sampleFreq);
    smoothKernel = 4;
    % COLLECT DATA IN THE PRESENCE OF VISUAL STIMULI
    Response = zeros(numChans,numStimuli,reps,stimLen);
    meanResponse = zeros(numChans,numStimuli,stimLen);
    for ii=1:numChans
        for jj=1:numStimuli
            stimStrobes = strobeTimes(svStrobed == jj);
            for kk=1:reps
                stimOnset = stimStrobes(kk);
                [~,index] = min(abs(timeStamps-stimOnset));
                temp = ChanData(index:index+stimLen-1,ii);
                Response(ii,jj,kk,:) = temp;
            end
            meanResponse(ii,jj,:) = smooth(mean(squeeze(Response(ii,jj,:,:)),1),smoothKernel);
        end
    end

end

function [significantStimuli,minLatency,maxLatency] = KStest(Response,meanResponse,baseStats)
        % KS TEST - distribution of VEP latencies significantly different 
        % in presence of a stimulus than a uniform distribution
        % use Benjamini-Hochberg correction for multiple comparisons
    global numChans numStimuli stimLen sampleFreq maxWin minWin;
   
    [~,minLatency] = min(Response(:,:,:,minWin),[],4);
    [~,maxLatency] = max(Response(:,:,:,maxWin),[],4);
    
    %baseWin = 1:round(0.02*sampleFreq);
    minVals = -min(meanResponse,[],3);
    significantStimuli = zeros(numChans,numStimuli);
    KSmins = zeros(numChans,numStimuli);
    KSmaxs = zeros(numChans,numStimuli);
    for ii=1:numChans
        for jj=1:numStimuli
            mintimes = squeeze(minLatency(ii,jj,:));
            maxtimes = squeeze(maxLatency(ii,jj,:));
            KSmins(ii,jj) = max(eCDF(mintimes,minWin)-linspace(0,1,length(minWin))');
            KSmaxs(ii,jj) = max(eCDF(maxtimes,maxWin)-linspace(0,1,length(maxWin))');
        end
        
        % add bootstrap step for confidence interval on KSmins and KSmaxs
        
%         [sortedPmins,~] = sort(squeeze(pmins(ii,:)));
%         l = (1:length(sortedPmins)).*alpha./length(sortedPmins);
%         threshLine = max(sortedPmins-l,0);
%         minSignif = find(threshLine==0,1,'last');
%         minthresh = min(sortedPmins(minSignif),alpha);
        

            for kk=1:numStimuli
%                 medianMin = median(minLatency(ii,kk,:));
%                 medianMax = median(maxLatency(ii,kk,:));
                if KSmins(ii,kk) > baseStats.minci(ii,2)
                    significantStimuli(ii,kk) = minVals(ii,kk);
                end
            end
    end
    minLatency = minLatency./sampleFreq;
    maxLatency = (maxLatency+maxWin(1)-1)./sampleFreq;
end

function [stat,se,ci] = Bootmin(Data,alpha,N,n)
%Bootmin.m
%   Bootstrap estimate of the standard error of a parameter or statistic, 
%    such as the sample mean or median.


if nargin < 2
    alpha = 0.05;
    N = 5000;
    n = length(Data);
elseif nargin < 3
    N = 5000;
    n = length(Data);
elseif nargin < 4
    n = length(Data);
end

Tboot = zeros(N,1);
ci = zeros(2,1);
for ii=1:N
    indeces = random('Discrete Uniform',n,[n,1]);
    temp = Data(indeces,:);
    [~,latency] = min(mean(temp,1));
    Tboot(ii) = latency;
end
stat = mean(Tboot);
se = std(Tboot);
ci(1) = quantile(Tboot,alpha/2);
ci(2) = quantile(Tboot,1-alpha/2);
end

function [stat,se,ci] = Bootmax(Data,alpha,N,n)
%Bootmin.m
%   Bootstrap estimate of the standard error of a parameter or statistic, 
%    such as the sample mean or median.
global maxWin;

if nargin < 2
    alpha = 0.05;
    N = 5000;
    n = length(Data);
elseif nargin < 3
    N = 5000;
    n = length(Data);
elseif nargin < 4
    n = length(Data);
end

Tboot = zeros(N,1);
ci = zeros(2,1);
for ii=1:N
    indeces = random('Discrete Uniform',n,[n,1]);
    temp = Data(indeces,maxWin);
    [~,latency] = max(mean(temp,1));
    Tboot(ii) = latency;
end
stat = mean(Tboot);
se = std(Tboot);
ci(1) = maxWin(1)+quantile(Tboot,alpha/2);
ci(2) = maxWin(1)+quantile(Tboot,1-alpha/2);
end

function [centerMass] = GetReceptiveField(significantStimuli,AnimalName)
    global numChans centerVals numStimuli;
    centerMass = struct('x',zeros(numChans,1),'y',zeros(numChans,1),...
        'Sigma',zeros(numChans,2,2));
    for ii=1:numChans
        dataX = [];dataY = [];
        for jj=1:numStimuli
            dataX = [dataX;repmat(centerVals(jj,1),[round(significantStimuli(ii,jj)),1])];
            dataY = [dataY;repmat(centerVals(jj,2),[round(significantStimuli(ii,jj)),1])];
            % this step is fairly odd. To fit the 2D Gaussian in units of pixels,
            % rather than in units of VEP magnitude, I make a distribution in
            % each dimension by creating a vector of pixel values. If the VEP
            % magnitude at the location (1200,300), (x,y), was 250 microVolts,
            % then the dataX vector will have 250 repetitions of the value 1200
            % and the dataY vector will have 250 repetitions of the value 300.
            % So, the distribution ends up being weighted more heavily by pixel
            % values with strong VEPs. The values obtained by this method for
            % the center of the retinotopic map are identical to those obtained
            % by performing a center of mass calculation (center of mass
            % exactly like those done in physics, SUM x*m / SUM m , except that
            % m is the VEP magnitude instead of the mass).
        end
        data = [dataX,dataY];

        try
            mnPDF = fitgmdist(data,1);
            centerMass.x(ii) = mnPDF.mu(1);
            centerMass.y(ii) = mnPDF.mu(2);
            centerMass.Sigma(ii,:,:) = mnPDF.Sigma;
        catch
            display(sprintf('\nError, Animal %d , Channel %d.\n',AnimalName,ii));
            centerMass.x(ii) = NaN;
            centerMass.y(ii) = NaN;
            centerMass.Sigma(ii,:,:) = NaN;
        end

    end
end

function [stimVals,h] = MakePlots(significantStimuli,meanResponse,centerMass,AnimalName,minInterval,maxInterval)
    global numChans numStimuli w_pixels h_pixels centerVals Radius stimLen;
    sigma = Radius;
    halfwidth = 3*sigma;
    [xx,yy] = meshgrid(-halfwidth:halfwidth,-halfwidth:halfwidth);
    gaussian = exp(-(xx.*xx+yy.*yy)./(2*sigma^2));
    stimVals = zeros(numChans,w_pixels,h_pixels);
    x=1:w_pixels;
    y=1:h_pixels;
    
    minVals = min(meanResponse,[],3);
    maxVals = max(meanResponse,[],3);
    xPos = unique(centerVals(:,1));
    yPos = unique(centerVals(:,2));
    
    xDiff = mean(diff(xPos));
    yDiff = mean(diff(yPos));
    xconv = stimLen/xDiff; % the max(diff(sort ...
                     % is equal to the width of the mapping stimulus (the width
                     % of the square that encloses the sinusoidal grating with
                     % overlain 2D Gaussian kernel) ... this value is
                     % equivalent to 2*Radius, which is output in the
                     % RetinoStim file, so I could use that. The only issue is
                     % that doing so would potentially preclude other types of
                     % stimuli from being analyzed by this code. As is, the
                     % code is fairly general to accept different types of
                     % mapping stimuli
    yconv = 1000/yDiff; % for height of the stimulus


    for ii=1:numChans
        h(ii) = figure;
    end

    for ii=1:numChans
        figure(h(ii));axis([0 w_pixels 0 h_pixels]);
        title(sprintf('VEP Retinotopy, Channel %d, Animal %d',ii,AnimalName));
        xlabel('Horizontal Screen Position (pixels)');ylabel('Vertical Screen Position (pixels)');
        hold on;

        for jj=1:numStimuli
            minErrors = round(minInterval(ii,jj,1)):round(minInterval(ii,jj,2));
            minErrLen = length(minErrors);
            maxErrors = round(maxInterval(ii,jj,1)):round(maxInterval(ii,jj,2));
            maxErrLen = length(maxErrors);
            tempx = centerVals(jj,1);
            tempy = centerVals(jj,2);
            stimVals(ii,tempx-Radius:tempx+Radius,tempy-Radius:tempy+Radius) = significantStimuli(ii,jj);
            plot(((1:1:stimLen)./xconv+centerVals(jj,1)-0.5*xDiff),...
                (squeeze(meanResponse(ii,jj,:))'./yconv+centerVals(jj,2)),'k','LineWidth',2);
            plot(minErrors./xconv+centerVals(jj,1)-0.5*xDiff,...
                (ones(minErrLen,1).*minVals(ii,jj))./yconv...
                +centerVals(jj,2)-10,'k','LineWidth',2);
            plot(maxErrors./xconv+centerVals(jj,1)-0.5*xDiff,...
                (ones(maxErrLen,1).*maxVals(ii,jj))./yconv...
                +centerVals(jj,2)+10,'k','LineWidth',2);
            plot((tempx-Radius)*ones(Radius*2+1,1),(tempy-Radius):(tempy+Radius),'k','LineWidth',2);
            plot((tempx+Radius)*ones(Radius*2+1,1),(tempy-Radius):(tempy+Radius),'k','LineWidth',2);
            plot((tempx-Radius):(tempx+Radius),(tempy-Radius)*ones(Radius*2+1,1),'k','LineWidth',2);
            plot((tempx-Radius):(tempx+Radius),(tempy+Radius)*ones(Radius*2+1,1),'k','LineWidth',2);

        end
        if isnan(centerMass.x(ii)) ~= 1
           temp = squeeze(stimVals(ii,:,:))';
%             blurred = conv2(temp,gaussian,'same');
%             blurred = (blurred./max(max(blurred))).*max(max(temp));
            imagesc(x,y,temp,'AlphaData',0.5);set(gca,'YDir','normal');w=colorbar;
            ylabel(w,'VEP Negativity (\muV)');colormap('jet');
            hold off;

        end
    end
end

function [Fx] = eCDF(Data,win)
%eCDF.m
%   Creation of the empirical distribution function (empirical cumulative
%   distribution function) for an array of Data
%
%INPUT: Data - the data as a vector
%       OPTIONAL:
%       alpha - confidence level for nonparametric 1-alpha confidence
%         bands, defaults to 0.05 for a 95% confidence band
%OUTPUT: Fx - the empirical distribution function
%        x - the points at which Fx is calculated

%Created: 2016/07/09
%  Byron Price
%Updated: 2016/09/07
%By: Byron Price

n = length(Data);
Fx = zeros(length(win),1);
Data = sort(Data);

count = 1;
for ii=win
    Fx(count) = sum(Data<=ii)/n;
    count = count+1;
end

end

