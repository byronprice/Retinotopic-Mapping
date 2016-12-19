function [] = RealTimeRetinoBayes_StimComp(AnimalName)
% RealTimeRetino_StimComp.m
%  Client/stimulus computer side code to perform a real-time, closed-loop 
%   retinopic mapping experiment. Will require a number of dependencies,
%   including code running on the recording computer
%   (RealTimeRetinoBayes_RecordComp.m), Psychtoolbox, and a function to
%   create a USB object to send TTL pulses from the stimulus computer to
%   the recording computer. To run, open Plexon, go to VEP settings for two
%   channels (if more than two, the code will work, but you'll need to
%   adjust for that), start data acquisition, run
%   RealTimeRetinoBayes_RecordComp on the recording computer in MATLAB,
%   then run this function.

% The algorithm performs Bayesian updating based on a likelihood function
%  calculated from hundreds of data points. That data gives the probability
%  of a VEP-like event at a specific location on the screen given the
%  distance from that location to the retinotopic center of mass. If the
%  location is itself the center of mass, then a VEP-like event occurs with
%  probability = 0.7 (approximately). If the location is not the center of
%  mass, then this value decays away exponentially to about 0.2 . Given
%  this knowledge, we can flash a stimulus and record either a HIT or a
%  MISS. Then, we can calculate the likelihood at each point in space given
%  the data of a HIT or a MISS. A HIT at one point in space tells us
%  something about that point, but it also tells us something about the
%  points nearby and far away. A MISS is similarly informative, for if we
%  miss for example in the top right corner 4/5 times, then we're pretty
%  certain the retinotopic center of mass is far away from there.

%Input: AnimalName - name of the animal, e.g. 12345

%Created: 2016/09/20, 24 Cummington Mall, Boston
% Byron Price
%Updated: 2016/12/19
%  By: Byron Price

cd('~/CloudStation/ByronExp/Retino');
load('BayesVars_Gauss.mat');

maxPrHit = b(1)+b(2);

% PARAMETERS for communication with recording computer
startEXP = 254;
endEXP = 255;

startRUN = 252;
endRUN = 253;


tcpipClient = tcpip('128.197.59.169',30000,'NetworkRole','client');
bufferSize = 1000; % bytes, (we won't need this much)
set(tcpipClient,'InputBufferSize',bufferSize);
set(tcpipClient,'Timeout',1);
fopen(tcpipClient);

directory = '~/Documents/MATLAB/Byron/Retinotopic-Mapping';
cd(directory);

global GL;

% Make sure this is running on OpenGL Psychtoolbox:
AssertOpenGL;

% usb = ttlInterfaceClass.getTTLInterface;
usb = usb1208FSPlusClass;
display(usb);

% Choose screen with maximum id - the secondary display:
screenid = max(Screen('Screens'));

% Open a fullscreen onscreen window on that display, choose a background
% color of 127 = gray with 50% max intensity; 0 = black; 255 = white
background = 127;
[win,~] = Screen('OpenWindow', screenid,background);

gammaTable = makeGrayscaleGammaTable(gama,0,255);
Screen('LoadNormalizedGammaTable',win,gammaTable);

% Switch color specification to use the 0.0 - 1.0 range
Screen('ColorRange', win, 1);

% Query window size in pixels
[w_pixels, h_pixels] = Screen('WindowSize', win);

% Retrieve monitor refresh duration
ifi = Screen('GetFlipInterval', win);

dgshader = [directory '/Retinotopy.vert.txt'];
GratingShader = LoadGLSLProgramFromFiles({ dgshader, [directory '/Retinotopy.frag.txt'] }, 1);
gratingTex = Screen('SetOpenGLTexture', win, [], 0, GL.TEXTURE_3D,w_pixels,...
    h_pixels, 1, GratingShader);

% screen size in millimeters and a conversion factor to get from mm to pixels
[w_mm,h_mm] = Screen('DisplaySize',screenid);
conv_factor = (w_mm/w_pixels+h_mm/h_pixels)/2;
mmPerPixel = conv_factor;
conv_factor = 1/conv_factor;

% perform unit conversions
Radius = (tan(degreeRadius*pi/180)*(DistToScreen*10))*conv_factor; % get number of pixels
     % that degreeRadius degrees of visual space will occupy
Radius = round(Radius);
temp = (tan((1/degreeSpatFreq)*pi/180)*(DistToScreen*10))*conv_factor;
spatFreq = 1/temp;

orientation = rand([repMax,1]).*pi;

xaxis = 1:10:w_pixels;
yaxis = 1:10:h_pixels;

numPositions = length(xaxis)*length(yaxis);
centerVals = zeros(numPositions,2);

% what follows for hitPriors and missPriors amounts to a compensation for
% the 2-D convolution of our likelihood functions, which will tend to
% accumulate probability mass in the center of the screen
stimSelection = ones(numPositions,1)./numPositions;

count = 1;
for ii=1:length(xaxis)
    for jj=1:length(yaxis)
        centerVals(count,1) = xaxis(ii);
        centerVals(count,2) = yaxis(jj);
%         if xaxis(ii) < 400
%             Prior(count,1) = 0;
%         elseif xaxis(ii) > w_pixels-400
%             Prior(count,2) = 0;
%         end
        count = count+1;
    end
end
DistFun = @(CenterPoint,AllPoints) (ceil(sqrt((CenterPoint(1)-AllPoints(:,1)).^2+(CenterPoint(2)-AllPoints(:,2)).^2))+1);

Response = zeros(numChans,repMax,2);

% Define first and second ring color as RGBA vector with normalized color
% component range between 0.0 and 1.0, based on Contrast between 0 and 1
% create all textures in the same window (win), each of the appropriate
% size
Grey = 0.5;
Black = 0;
White = 1;

Screen('BlendFunction',win,GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);

stimTimes = stimTime-0.1+exprnd(0.1,[repMax,1]);
WaitSecs(startPause);
usb.strobeEventWord(startEXP);
PrHitNoise = fread(tcpipClient,numChans,'double');

Priority(9);
% Mapping Loop
vbl = Screen('Flip',win);

PrevIndex = 1;
count = 1;
while count < repMax %&& max(Prior(:,ii)) < thresh)
    dists = DistFun(centerVals(PrevIndex,:),centerVals);
    unifRand = rand;
    tempDist = stimSelection(:,1);
    tempDist(dists<Radius) = 0;
    tempDist = tempDist./sum(tempDist);
    CDF = cumsum(tempDist(:));
    CDF = CDF-min(CDF);
    temp = unifRand-CDF;
    temp(temp<0) = 0;
    [~,index] = min(temp);
    stimCenter = centerVals(index,:);
    PrevIndex = index;
    vbl = Screen('Flip',win);
    usb.strobeEventWord(startRUN);
    
    % Draw the procedural texture as any other texture via 'DrawTexture'
    Screen('DrawTexture', win,gratingTex, [],[],...
        [],[],[],[Grey Grey Grey Grey],...
        [], [],[White,Black,...
        Radius,stimCenter(1),stimCenter(2),spatFreq,orientation(count),0]);
    % Request stimulus onset
    vbl = Screen('Flip', win,vbl+ifi/2);usb.strobeEventWord(1);
    vbl = Screen('Flip',win,vbl-ifi/2+stimTimes(count));
    
    WaitSecs(0.2);
    usb.strobeEventWord(endRUN);
    
    WaitSecs(0.3);
    data = fread(tcpipClient,numChans,'double'); % hit or miss on each channel
    if isempty(data) == 0
            Response(:,count,1) = index*ones(numChans,1);
            Response(:,count,2) = data;

            % Bayesian update step
            
            %             if hit_or_miss == 1
            %                 Posterior = Likelihood_Hit(allPossDists)'.*Prior(:,ii);
            %                 Prior(:,ii) = Posterior./sum(Posterior);
            %                 display('HIT');
            %             elseif hit_or_miss == 0
            %                 Posterior = Likelihood_Miss(allPossDists)'.*Prior(:,ii);
            %                 Prior(:,ii) = Posterior./sum(Posterior);
            %                 display('MISS');
            %             end
    else
        continue;
    end
    count = count+1;
    clear data;
end
usb.strobeEventWord(endEXP);

% remove values that failed to register
Response = Response(:,squeeze(Response(1,:,1))>0,:);
reps = size(Response,2);

% calculate likelihood, which is a model with two parameters, retinotopic 
%  center of mass and standard deviation, sigma
maxDist = ceil(sqrt((xaxis(end)-xaxis(1)).^2+(yaxis(end)-yaxis(1)).^2));
DistToCenterMass = (0:maxDist)';

maxSigma = 800;
Sigma = 50:maxSigma;numSigmas = length(Sigma);

numDists = length(DistToCenterMass);
SigmaSquare = 1./(Sigma.^2);
% DistToCenterMass must be a column vector, SigmaSquare must be a row vector
Likelihood = zeros(numDists,numSigmas,numChans);
for ii=1:numChans
    tempB = b;
    tempB(2) = PrHitNoise(ii);
    tempB(1) = min(maxPrHit,tempB(2)+b(1));
    Likelihood(:,:,ii) = likelihoodFun(tempB,DistToCenterMass,SigmaSquare);
end

TotalLikelihood = ones(numPositions,numSigmas,numChans);
%allPossDists = zeros(numPositions,1);
stimCenter = zeros(2,1);
for ii=1:reps
    index = Response(1,ii,1);
    stimCenter(1) = centerVals(index,1);stimCenter(2) = centerVals(index,2);
    
    allPossDists = DistFun(stimCenter,centerVals);
    for kk=1:numChans
        hit_or_miss = Response(kk,ii,2);
        TotalLikelihood(:,:,kk) = squeeze(TotalLikelihood(:,:,kk)).*(squeeze(Likelihood(allPossDists,:,kk)).^hit_or_miss)...
            .*((1-squeeze(Likelihood(allPossDists,:,kk))).^(1-hit_or_miss));
    end
end
      
% massive Posterior matrix, number of center positions by number of
%  possible values for the standard deviation sigma
Posterior = zeros(numPositions,numSigmas,numChans);

h(1) = figure();
h(2) = figure();
h(3) = figure();
stimVals = zeros(numChans,length(xaxis),length(yaxis));
for ii=1:numChans
    Posterior(:,:,ii) = TotalLikelihood(:,:,ii); % multiplied by two priors based on 
                        % the above inverse convolutions 
    Posterior(:,:,ii) = Posterior(:,:,ii)./sum(sum(Posterior(:,:,ii)));
    
    marginalPost = sum(squeeze(Posterior(:,:,ii)),2);
    count = 1;
    for jj=1:length(xaxis)
        for kk=1:length(yaxis)
            stimVals(ii,jj,kk) = marginalPost(count);
            count = count+1;
        end
    end
    figure(h(1));subplot(numChans,1,ii);
    imagesc(xaxis,yaxis,squeeze(stimVals(ii,:,:)./(1/numPositions))');
    set(gca,'YDir','normal');caxis([0,50]);
    w = colorbar;ylabel(w,'Normalized Probability (1 = Pr under Discrete Uniform)');
    title(sprintf('VEP Retinotopy [Pr(Given Point is Center of Mass)]: Channel %d, Animal %d',ii,AnimalName));
    xlabel('Horizontal Screen Position');ylabel('Vertical Screen Position');
    
    figure(h(2));subplot(numChans,1,ii);
    imagesc(xaxis,yaxis,log(squeeze(stimVals(ii,:,:))'));
    set(gca,'YDir','normal');caxis([-20 -1]);
    w = colorbar;ylabel(w,'Natural Log Probability');
    title(sprintf('VEP Retinotopy [ln{Pr(Given Point is Center of Mass)}]: Channel %d, Animal %d',ii,AnimalName));
    xlabel('Horizontal Screen Position');ylabel('Vertical Screen Position');
    
    figure(h(3));subplot(numChans,1,ii);
    inds = squeeze(Response(ii,:,1));
    hitsMisses = squeeze(Response(ii,:,2));
    xVals = centerVals(inds,1);
    yVals = centerVals(inds,2);
    
    scatter(xVals(hitsMisses==1),yVals(hitsMisses==1),'ob');
    axis([0 max(xaxis) 0 max(yaxis)]);hold on;
    scatter(xVals(hitsMisses==0),yVals(hitsMisses==0),'xr');hold off;
end

Priority(0);
Screen('CloseAll');
cd('~/CloudStation/ByronExp/Retino');
Date = datetime('today','Format','yyyy-MM-dd');
Date = char(Date); Date = strrep(Date,'-','');Date=str2double(Date);

% save data
save(sprintf('BayesRetinoMap%d_%d.mat',Date,AnimalName),'centerVals',...
    'Posterior','mmPerPixel','degreeRadius','Response','degreeSpatFreq','numChans',...
    'repMax','stimTime','PrHitNoise','xaxis','yaxis','reps','numPositions','b');
end

function gammaTable = makeGrayscaleGammaTable(gamma,blackSetPoint,whiteSetPoint)
% Generates a 256x3 gamma lookup table suitable for use with the
% psychtoolbox Screen('LoadNormalizedGammaTable',win,gammaTable) command
% 
% gammaTable = makeGrayscaleGammaTable(gamma,blackSetPoint,whiteSetPoint)
%
%   gamma defines the level of gamma correction (1.8 or 2.2 common)
%   blackSetPoint should be the highest value that results in a non-unique
%   luminance value on the monitor being used (sometimes values 0,1,2, all
%   produce the same black pixel value; set to zero if this is not a
%   concern)
%   whiteSetPoint should be the lowest value that returns a non-unique
%   luminance value (deal with any saturation at the high end)
% 
%   Both black and white set points should be defined on a 0:255 scale

gamma = max([gamma 1e-4]); % handle zero gamma case
gammaVals = linspace(blackSetPoint/255,whiteSetPoint/255,256).^(1./gamma);
gammaTable = repmat(gammaVals(:),1,3);
end
