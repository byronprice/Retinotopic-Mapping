function [] = Retinotopy(AnimalName,Hemisphere,DistToScreen,degreeRadius)
%Retinotopy.m
%  Display a series of flashing squares to determine retinotopy of
%   LFP recording electrode.
%  Each square will occupy a 7-degree radius of visual space
% INPUT: Obligatory-
%        AnimalName - animal's unique identifier as a number, e.g. 45602
%
%        Optional- 
%        Hemisphere - hemisphere where electrodes are placed, 'LH' for
%           left, 'RH' for right, 'both' for both, defaults to 'both'
%        DistToScreen - physical distance of observer from the screen, in
%           units of cm
%        degreeRadius - degrees of visual field that radius of square will occupy
%        spatFreq - spatial frequency of oriented sinusoidal grating
%
% OUTPUT: a file with stimulus parameters named RetinoStimDate_AnimalName
%           e.g. RetinoStim20160708_12345.mat to be saved in the RetinoExp
%           folder under '/MATLAB/Byron/'
% Created: 2016/05/24 at 24 Cummington, Boston
%  Byron Price
% Updated: 2016/07/22
%  By: Byron Price

directory = '/home/jglab/Documents/MATLAB/Byron/Retinotopic-Mapping/';
if nargin < 2
    Hemisphere = 'both';
    DistToScreen = 25;
    degreeRadius = 5;
    reps = 40;
    stimLen = 50/1000;
    waitTime = 0.5;
    startPause = 120; % 120 seconds of silence before commencing
    spatFreq = 0.3;
elseif nargin < 3
    DistToScreen = 25;
    degreeRadius = 5;
    reps = 40;
    stimLen = 50/1000;
    waitTime = 0.5;
    startPause = 120; % 120 seconds of silence before commencing
    spatFreq = 0.3;
end

Date = datetime('today','Format','yyyy-MM-dd');
Date = char(Date); Date = strrep(Date,'-','');
% Acquire a handle to OpenGL, so we can use OpenGL commands in our code:
global GL;

% Make sure this is running on OpenGL Psychtoolbox:
AssertOpenGL;

usb = usb1208FSPlusClass;
display(usb);

WaitSecs(10);

% Choose screen with maximum id - the secondary display:
screenid = max(Screen('Screens'));

% Open a fullscreen onscreen window on that display, choose a background
% color of 128 = gray with 50% max intensity; 0 = black
[win,~] = Screen('OpenWindow', screenid,0);

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
temp = (tan((1/spatFreq)*pi/180)*(DistToScreen*10))*conv_factor;
spatFreq = 1/temp;

if strcmp(Hemisphere,'LH') == 1
    centerX = round(w_pixels/2)-100:2*Radius:w_pixels;
    centerY = Radius+1:2*Radius:h_pixels-Radius/2;
elseif strcmp(Hemisphere,'RH') == 1
    centerX = Radius+1:2*Radius:round(w_pixels/2)+100;
    centerY = Radius+1:2*Radius:h_pixels-Radius/2;
elseif strcmp(Hemisphere,'both') == 1
    centerX = 3*Radius:2*Radius:w_pixels-3*Radius;
    centerY = Radius+1:2*Radius:h_pixels-Radius/2;
end
numStimuli = length(centerX)*length(centerY);

centerVals = zeros(numStimuli,2);
count = 1;
for ii=1:length(centerX)
    for jj=1:length(centerY)
        centerVals(count,1) = centerX(ii);
        centerVals(count,2) = centerY(jj);
        count = count+1;
    end
end

for ii=1:50
    indeces = randperm(numStimuli,numStimuli);
    centerVals = centerVals(indeces,:);
end

estimatedTime = ((waitTime+stimLen)*reps*numStimuli+startPause)/60;
display(sprintf('Estimated time: %3.2f minutes',estimatedTime));


% Define first and second ring color as RGBA vector with normalized color
% component range between 0.0 and 1.0, based on Contrast between 0 and 1
% create all textures in the same window (win), each of the appropriate
% size
Grey = 0.5;
Black = 0;
White = 1;

%orientation = (rand([numStimuli,1])*360)*pi/180;
% Perform initial flip to gray background and sync us to the retrace:
Priority(9);

usb.startRecording;
WaitSecs(startPause);

% Animation loop
for zz = 1:reps
    vbl = Screen('Flip', win);
    for ii=1:numStimuli
        orient = rand*2*pi;
        % Draw the procedural texture as any other texture via 'DrawTexture'
        Screen('DrawTexture', win,gratingTex, [],[],...
            [],[],[],[Grey Grey Grey Grey],...
            [], [],[White,Black,...
            Radius,centerVals(ii,1),centerVals(ii,2),spatFreq,orient,0]);
        % Request stimulus onset
        vbl = Screen('Flip', win, vbl + ifi/2);usb.strobe;
        vbl = Screen('Flip',win,vbl+ifi/2+stimLen);
        vbl = Screen('Flip',win,vbl+ifi/2+waitTime);
    end
    WaitSecs(2);
end
WaitSecs(2);
usb.stopRecording;
Priority(0);

cd('~/Documents/MATLAB/Byron/RetinoExp')
fileName = strcat('RetinoStim',Date,'_',num2str(AnimalName),'.mat');
save(fileName,'centerVals','Radius','reps','stimLen','startPause',...
    'numStimuli','w_pixels','h_pixels','spatFreq','mmPerPixel')
% Close window
Screen('CloseAll');

end
