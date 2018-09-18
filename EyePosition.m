function [] = EyePosition(AnimalName)
%EyePosition.m
%  Display a series of flashing gabors to determine where the mouse is
%   looking
%  Each circle will occupy an ~ 2.5-degree radius of visual space
% INPUT: Obligatory-
%        AnimalName - animal's unique identifier as a number, e.g. 45602
%  
%
% OUTPUT: a file with stimulus parameters named EyePosStimDate_AnimalName
%           e.g. EyePosStim20160708_12345.mat to be saved in the Retino
%           folder on the CloudStation
% Created: 2018/08/22 at 24 Cummington, Boston
%  Byron Price
% Updated: 2018/08/22
%  By: Byron Price

currentdirectory = '~/Documents/MATLAB/Byron/Retinotopic-Mapping';
cd(currentdirectory);

blocks = 20;
directions = 8;
numStimuli = 4;
holdTime = 30;
gama = 2.1806;
degreeRadius = 2.5;
DistToScreen = 25;
spatFreq = 0.01;
stimTime = 0.05;
delayTime = 0.2;
waitTime = 1;
stimDist = 10;

Date = datetime('today','Format','yyyy-MM-dd');
Date = char(Date); Date = strrep(Date,'-','');Date=str2double(Date);
% Acquire a handle to OpenGL, so we can use OpenGL commands in our code:
global GL;

% Make sure this is running on OpenGL Psychtoolbox:
AssertOpenGL;

% usb = ttlInterfaceClass.getTTLInterface;
usb = usb1208FSPlusClass;
display(usb);

WaitSecs(1);

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

dgshader = [currentdirectory '/Retinotopy.vert.txt'];
GratingShader = LoadGLSLProgramFromFiles({ dgshader, [currentdirectory '/Retinotopy.frag.txt'] }, 1);
gratingTex = Screen('SetOpenGLTexture', win, [], 0, GL.TEXTURE_3D,w_pixels,...
    h_pixels, 1, GratingShader);

% screen size in millimeters and a conversion factor to get from mm to pixels
[w_mm,h_mm] = Screen('DisplaySize',screenid);
conv_factor = (w_mm/w_pixels+h_mm/h_pixels)/2;
mmPerPixel = conv_factor;
conv_factor = 1/conv_factor;

% perform unit conversions
Radius = (tan((degreeRadius/2)*pi/180)*(DistToScreen*10*2))*conv_factor; % get number of pixels
     % that degreeRadius degrees of visual space will occupy
     
temp = (tan(((1/spatFreq)/2)*pi/180)*(DistToScreen*10*2))*conv_factor;
newSpatFreq = 1/temp;

stimDistPix = (tan((stimDist/2)*pi/180)*(DistToScreen*10*2))*conv_factor; % convert to pixels

centerPos = [w_pixels/2,h_pixels/2];
centerVals = zeros(directions,numStimuli,2);
stimStrobeNums = zeros(directions,numStimuli);

initialAngle = 0;
count = 2;
for ii=1:directions
    for jj=1:numStimuli
        stimStrobeNums(ii,jj) = count;
        centerVals(ii,jj,1) = centerPos(1)+(jj-1)*stimDistPix*cos(initialAngle*pi/180);
        centerVals(ii,jj,2) = centerPos(2)+(jj-1)*stimDistPix*sin(initialAngle*pi/180);
        count = count+1;
    end
    initialAngle = initialAngle+360/directions;
end

% estimatedTime = ((waitTime+0.05+stimTime)*reps*blocks+blocks*holdTime)/60;
fprintf('\nEstimated time: %3.2f minutes\n',13.3);

% Define first and second ring color as RGBA vector with normalized color
% component range between 0.0 and 1.0, based on Contrast between 0 and 1
% create all textures in the same window (win), each of the appropriate
% size
Grey = 0.5;
Black = 0;
White = 1;

phase = binornd(1,0.5,[blocks,directions]);
phase = 2.*phase-1;
phase = phase.*(pi/3);

Screen('BlendFunction',win,GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);

orientation = rand([blocks,directions]).*(2*pi);
waitTimes = waitTime+(rand([blocks,directions])*0.5-0.25);

% Perform initial flip to gray background and sync us to the retrace:
Priority(9);

usb.startRecording;WaitSecs(1);usb.strobeEventWord(1);
WaitSecs(holdTime*2);
% usb.strobeEventWord(startEXP);WaitSecs(1);

% Animation loop
for yy = 1:blocks
    for zz=1:directions
        Screen('DrawTexture',win,gratingTex, [],[],...
            [],[],[],[Grey Grey Grey Grey],...
            [], [],[White,Black,...
            Radius,centerVals(zz,1,1),centerVals(zz,1,2),newSpatFreq,orientation(yy,zz),phase(yy,zz)]);
        % Request stimulus onset
        vbl = Screen('Flip', win);
        usb.strobeEventWord(stimStrobeNums(zz,1));
        vbl = Screen('Flip',win,vbl-ifi/2+stimTime);
        vbl = Screen('Flip',win,vbl-ifi/2+delayTime*5);
        ii=1;
        while ii<=numStimuli
            % Draw the procedural texture as any other texture via 'DrawTexture'
            Screen('DrawTexture',win,gratingTex, [],[],...
                [],[],[],[Grey Grey Grey Grey],...
                [], [],[White,Black,...
                Radius,centerVals(zz,ii,1),centerVals(zz,ii,2),newSpatFreq,orientation(yy,zz),phase(yy,zz)]);
            % Request stimulus onset
            vbl = Screen('Flip', win);
            usb.strobeEventWord(stimStrobeNums(zz,ii));
            vbl = Screen('Flip',win,vbl-ifi/2+stimTime);
            vbl = Screen('Flip',win,vbl-ifi/2+delayTime);
            ii=ii+1;
        end
        vbl = Screen('Flip',win,vbl-ifi/2+waitTimes(yy,zz));
    end
    usb.strobeEventWord(1);
    vbl = Screen('Flip',win,vbl-ifi/2+holdTime/3);
end
% WaitSecs(1);usb.strobeEventWord(endEXP);
WaitSecs(holdTime);
usb.stopRecording;
Priority(0);

cd('~/CloudStation/ByronExp/Retino');

fileName = sprintf('EyePosStim%d_%d.mat',Date,AnimalName);
save(fileName,'centerVals','Radius','degreeRadius','stimStrobeNums','stimTime',...
    'holdTime','waitTimes','delayTime','blocks','directions','numStimuli',...
    'w_pixels','h_pixels','spatFreq','mmPerPixel','DistToScreen','orientation',...
    'phase')
% Close window
Screen('CloseAll');

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
