function [] = RealTimeRetino_StimComp(AnimalName)
% RealTimeRetino_StimComp.m
%  Client/stimulus computer side code to perform a real-time, closed-loop 
%   retinopic mapping experiment. Will require a number of dependencies,
%   including code running on the recording computer
%   (RealTimeRetino_RecordingComp.m), Psychtoolbox, and a function to
%   create a USB object to send TTL pulses from the stimulus computer to
%   the recording computer.

%Input: AnimalName - name of the animal, e.g. 12345

%Created: 2016/09/13, 24 Cummington Mall, Boston
% Byron Price
%Updated: 2016/09/15
%  By: Byron Price

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
set(tcpipClient,'Timeout',2);
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

degreeRadii = [25,20,15,10,5];
numTests = length(degreeRadii);
degreeSpatFreq = 0.05;

% perform unit conversions
Radii = (tan(degreeRadii*pi/180)*(DistToScreen*10))*conv_factor; % get number of pixels
     % that degreeRadius degrees of visual space will occupy
temp = (tan((1/degreeSpatFreq)*pi/180)*(DistToScreen*10))*conv_factor;
spatFreq = 1/temp;
orient = pi/4;

numStimuli = 4; % must be a perfect square
xaxis = w_pixels/(numStimuli-1);
yaxis = h_pixels/(numStimuli-1);
x = 1:w_pixels;y=1:h_pixels;
stimVals = zeros(numChans,w_pixels,h_pixels);

stimTime = 0.15;
WaitTime = 0.15;
repMax = 20;
numFlashes = 3;

strobeValues = 1:numTests*numStimuli;
% Define first and second ring color as RGBA vector with normalized color
% component range between 0.0 and 1.0, based on Contrast between 0 and 1
% create all textures in the same window (win), each of the appropriate
% size
Grey = 0.5;
Black = 0;
White = 1;

Screen('BlendFunction',win,GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);

WaitSecs(20);

usb.strobeEventWord(startEXP);
Pr = fread(tcpipClient,numChans,'double');

binoThresh = zeros(numChans,repMax);
binoThresh(:,1:3) = 3*numFlashes+1;
alpha = 0.05;
for ii=1:numChans
    for jj=4:repMax
        x = 1:jj*numFlashes;
        y = binopdf(x,jj*numFlashes,Pr(ii));
        [~,ind] = max(y);
        y(1:ind) = 1;
        Thresh = find(y<alpha,1,'first');
        binoThresh(ii,jj) = Thresh;
    end
end
binoThresh

Priority(9);
% Mapping Loop
vbl = Screen('Flip',win);
for ii=1:numChans
    Data = zeros(numTests*numStimuli,1);
    centerVals = zeros(numStimuli,2);
    
    ccount = 1;
    for ii=1:sqrt(numStimuli)
        for jj=1:sqrt(numStimuli)
            centerVals(ccount,1) = xaxis*ii;
            centerVals(ccount,2) = yaxis*jj;
            ccount = ccount+1;
        end
    end
    for jj=1:numTests
        check = 0;
        count = 1;
        while check == 0
            vbl = Screen('Flip',win);
            indeces = randperm(numStimuli);
            usb.strobeEventWord(startRUN);
            for ll=1:numStimuli
                for nn=1:numFlashes
                    % Draw the procedural texture as any other texture via 'DrawTexture'
                    Screen('DrawTexture', win,gratingTex, [],[],...
                        [],[],[],[Grey Grey Grey Grey],...
                        [], [],[White,Black,...
                        Radii(jj),centerVals(indeces(ll),1),centerVals(indeces(ll),2),spatFreq,orient,0]);
                    % Request stimulus onset
                    vbl = Screen('Flip', win,vbl+ifi/2);usb.strobeEventWord(strobeValues((jj-1)*numStimuli+indeces(ll)));
                    vbl = Screen('Flip',win,vbl-ifi/2+stimTime);
                    vbl = Screen('Flip',win,vbl-ifi/2+WaitTime);
                end
                vbl = Screen('Flip',win,vbl-ifi/2+stimTime*2);
            end
            usb.strobeEventWord(endRUN);
            dataSize = fread(tcpipClient,3,'double');
            
            if isempty(dataSize) == 0
                temp = fread(tcpipClient,dataSize(1)*dataSize(2),'double');
                temp = reshape(temp,[dataSize(1),dataSize(2)]);
                for mm=1:dataSize(1)
                    index = temp(mm,2);
                    if index <= numel(Data)
                        Data(index) = Data(index)+temp(mm,1);
                        newInd = index-(jj-1)*numStimuli;
                        lowx = max(round(centerVals(newInd,1)-Radii(jj)),1);
                        highx = min(round(centerVals(newInd,1)+Radii(jj)),w_pixels);
                        lowy = max(round(centerVals(newInd,2)-Radii(jj)),1);
                        highy = min(round(centerVals(newInd,2)+Radii(jj)),h_pixels);
                        stimVals(ii,lowx:highx,lowy:highy) = stimVals(ii,lowx:highx,lowy:highy)+1;
                    end
                end
            end
            check = sum(Data((jj-1)*numStimuli+1:end) >= binoThresh(ii,count));
            count = count+1;
            if count == repMax
                check = 1;
            end
        end

        if jj < numTests
            [~,winIndex] = max(Data((jj-1)*numStimuli+1:end));
            baseCenter = centerVals(winIndex,:);
            for mm=1:numStimuli
                centerVals(mm,1) = max(baseCenter(1)+cos((pi/4)*((mm-1)*2+1))*Radii(jj+1),1);
                centerVals(mm,2) = max(baseCenter(2)+sin((pi/4)*((mm-1)*2+1))*Radii(jj+1),1);
            end
        end
    end
    usb.strobeEventWord(endCHAN);
end
figure();imagesc(x,y,stimVals(1,:,:)');

Priority(0);
Screen('CloseAll');
cd('~/CloudStation/ByronExp/Retino');
Date = datetime('today','Format','yyyy-MM-dd');
Date = char(Date); Date = strrep(Date,'-','');Date=str2double(Date);

% save data
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