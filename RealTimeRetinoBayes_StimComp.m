function [] = RealTimeRetinoBayes_StimComp(AnimalName)
% RealTimeRetino_StimComp.m
%  Client/stimulus computer side code to perform a real-time, closed-loop 
%   retinopic mapping experiment. Will require a number of dependencies,
%   including code running on the recording computer
%   (RealTimeRetino_RecordingComp.m), Psychtoolbox, and a function to
%   create a USB object to send TTL pulses from the stimulus computer to
%   the recording computer.

%Input: AnimalName - name of the animal, e.g. 12345

%Created: 2016/09/17, 24 Cummington Mall, Boston
% Byron Price
%Updated: 2016/09/17
%  By: Byron Price

cd('~/CloudStation/ByronExp/Retino');
load('BayesVars.mat');


% INITIAL PARAMETERS
startEXP = 254;
endEXP = 255;

startRUN = 252;
endRUN = 253;

endCHAN = 251;

numChans = 2;

tcpipClient = tcpip('128.197.59.169',30000,'NetworkRole','client');
bufferSize = 50000; % bytes, a big number (won't need this much)
set(tcpipClient,'InputBufferSize',bufferSize);
set(tcpipClient,'Timeout',1);
fopen(tcpipClient);

directory = '~/Documents/MATLAB/Byron/Retinotopic-Mapping';

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

gama = 2.1806;
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

DistToScreen = 25;

degreeRadius = 10;
degreeSpatFreq = 0.05;

% perform unit conversions
Radius = (tan(degreeRadius*pi/180)*(DistToScreen*10))*conv_factor; % get number of pixels
     % that degreeRadius degrees of visual space will occupy
Radius = round(Radius);
temp = (tan((1/degreeSpatFreq)*pi/180)*(DistToScreen*10))*conv_factor;
spatFreq = 1/temp;
orient = pi/6;

xaxis = Radius:10:w_pixels-Radius;
yaxis = 1:10:h_pixels-Radius;

numStimuli = length(xaxis)*length(yaxis);
centerVals = zeros(numStimuli,2);

Prior = ones(numStimuli,numChans)./numStimuli;
thresh = min(Prior(1,1)*100,0.05);

count = 1;
for ii=1:length(xaxis)
    for jj=1:length(yaxis)
        centerVals(count,1) = xaxis(ii);
        centerVals(count,2) = yaxis(jj);
        count = count+1;
    end
end

maxDist = ceil(sqrt((xaxis(end)-xaxis(1)).^2+(yaxis(end)-yaxis(1)).^2));
DistToCenterMass = 0:maxDist;

Likelihood_Hit = likelihoodFun(b,DistToCenterMass);
Likelihood_Miss = 1-Likelihood_Hit;

allPossDists = zeros(numStimuli,1);

stimTime = 0.05;
repMax = 100;

% Define first and second ring color as RGBA vector with normalized color
% component range between 0.0 and 1.0, based on Contrast between 0 and 1
% create all textures in the same window (win), each of the appropriate
% size
Grey = 0.5;
Black = 0;
White = 1;

Screen('BlendFunction',win,GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);

WaitSecs(15);

usb.strobeEventWord(startEXP);
Pr = fread(tcpipClient,numChans,'double');

Priority(9);
% Mapping Loop
vbl = Screen('Flip',win);
for ii=1:numChans  
    count = 1;
    while (count < repMax && max(Prior(:,ii)) < thresh)
        unifRand = rand;
        CDF = cumsum(Prior(:,ii));
        CDF = CDF-min(CDF);
        temp = unifRand-CDF;
        temp(temp<0) = 0;
        [~,index] = min(temp);
        stimCenter = centerVals(index,:);
        vbl = Screen('Flip',win);
        usb.strobeEventWord(startRUN);

        % Draw the procedural texture as any other texture via 'DrawTexture'
        Screen('DrawTexture', win,gratingTex, [],[],...
            [],[],[],[Grey Grey Grey Grey],...
            [], [],[White,Black,...
            Radius,stimCenter(1),stimCenter(2),spatFreq,orient,0]);
        % Request stimulus onset
        vbl = Screen('Flip', win,vbl+ifi/2);usb.strobeEventWord(1);
        vbl = Screen('Flip',win,vbl-ifi/2+stimTime);
        WaitSecs(abs(0.2-stimTime));
        usb.strobeEventWord(endRUN);
        
        data = fread(tcpipClient,2,'double');
        if isempty(data) == 0
            hit_or_miss = data(1);
            for jj=1:numStimuli
                allPossDists(jj) = ceil(sqrt((stimCenter(1)-centerVals(jj,1)).^2+...
                    (stimCenter(2)-centerVals(jj,2)).^2))+1;
            end
            % Bayesian update step
            if hit_or_miss == 1
                Posterior = Likelihood_Hit(allPossDists)'.*Prior(:,ii);
                Prior(:,ii) = Posterior./sum(Posterior);
            elseif hit_or_miss == 0
                Posterior = Likelihood_Miss(allPossDists)'.*Prior(:,ii);
                Prior(:,ii) = Posterior./sum(Posterior);
            end
        else
            continue;
        end

        count = count+1;
        clear data;
    end
    usb.strobeEventWord(endCHAN);
end

Posterior = Prior;
Priority(0);
Screen('CloseAll');
cd('~/CloudStation/ByronExp/Retino');
Date = datetime('today','Format','yyyy-MM-dd');
Date = char(Date); Date = strrep(Date,'-','');Date=str2double(Date);

% save data
save(sprintf('BayesRetinoMap%d_%d.mat',Date,AnimalName),'centerVals','Posterior');
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