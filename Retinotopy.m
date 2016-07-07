function [] = Retinotopy(AnimalName,Hemisphere,DistToScreen,degreeRadius)
%Retinotopy.m
%  Display a series of flashing circles to determine retinotopy of
%   LFP recording electrode.
%  Each circle will occupy a 4-degree radius of visual space
% INPUT: DistToScreen - physical distance of observer from the screen, in
%           units of cm
%        AnimalName - animal's unique identifier as a number, e.g. 45602
%        degreeRadius - degrees of visual field that radius of circle will occupy
%        Hemisphere - hemisphere where electrodes are placed, 'LH' for
%           left, 'RH' for right, 'both' for both, defaults to 'LH'
%
% Created: 2016/05/24 at 24 Cummington, Boston
%  Byron Price
% Updated: 2016/07/01
%  By: Byron Price

directory = pwd;
if nargin < 2
    Hemisphere = 'LH';
    DistToScreen = 25;
    degreeRadius = 5;
    checkSize = 5;
    reps = 20;
    stimLen = 50/1000;
    waitTime = 1;
    startPause = 60; % 60 seconds of silence before commencing
    numStimuli = 100;
elseif nargin < 3
    DistToScreen = 25;
    degreeRadius = 5;
    reps = 20;
    stimLen = 50/1000;
    waitTime = 1;
    startPause = 60; % 60 seconds of silence before commencing
    numStimuli = 100;
end

Date = datetime('today','Format','yyyy-MM-dd');
Date = char(Date); Date = strrep(Date,'-','');
% Acquire a handle to OpenGL, so we can use OpenGL commands in our code:
global GL;

% Make sure this is running on OpenGL Psychtoolbox:
AssertOpenGL;

usb = usb1208FSPlusClass;
display(usb);
WaitSecs(5);

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

% screen size in millimeters and a conversion factor to get from mm to pixels
[w_mm,h_mm] = Screen('DisplaySize',screenid);
conv_factor = (w_mm/w_pixels+h_mm/h_pixels)/2;

% perform unit conversions
Radius = (tan(degreeRadius*(2*pi)/360)*(DistToScreen*10))./conv_factor; % get number of pixels
     % that degreeRadius degrees of visual space will occupy
Radius = round(Radius);
checkSize = (tan(checkSize*(2*pi)/360)*(DistToScreen*10))./conv_factor;
checkSize = round(checkSize);

if strcmp(Hemisphere,'LH') == 1
    centerX = round(w_pixels/2)-50:2*Radius:w_pixels-Radius;
    centerY = Radius+1:2*Radius:h_pixels-Radius;
elseif strcmp(Hemisphere,'RH') == 1
    centerX = Radius+1:2*Radius:round(w_pixels/2)+50;
    centerY = Radius+1:2*Radius:h_pixels-Radius;
elseif strcmp(Hemisphere,'both') == 1
    centerX = Radius+1:2*Radius:w_pixels-Radius;
    centerY = Radius+1:2*Radius:h_pixels-Radius;
end

centerVals = zeros(numStimuli,2);
for ii=1:numStimuli
    centerVals(ii,1) = centerX(randi([1,numel(centerX)],1));
    centerVals(ii,2) = centerY(randi([1,numel(centerY)],1));
end
%permIndeces = randperm(length(centerVals));
%centerVals = centerVals(permIndeces,:);

dgshader = [directory '/Retinotopy.vert.txt'];
GratingShader = LoadGLSLProgramFromFiles({ dgshader, [directory '/Retinotopy.frag.txt'] }, 1);
gratingTex = Screen('SetOpenGLTexture', win, [], 0, GL.TEXTURE_3D,w_pixels,...
    h_pixels, 1, GratingShader);

% Define first and second ring color as RGBA vector with normalized color
% component range between 0.0 and 1.0, based on Contrast between 0 and 1
% create all textures in the same window (win), each of the appropriate
% size
Grey = 0.5;
Black = 0;
White = 1;

% Perform initial flip to gray background and sync us to the retrace:
vbl = Screen('Flip', win);

usb.startRecording;
WaitSecs(startPause);

% Animation loop
for zz = 1:reps
    for ii=1:length(centerVals)
        random = round(rand*checkSize);
        % Draw the procedural texture as any other texture via 'DrawTexture'
        Screen('DrawTexture', win,gratingTex, [],[],...
            [],[],[],[Grey Grey Grey Grey],...
            [], [],[White,Black,...
            Radius,centerVals(ii,1),centerVals(ii,2),checkSize,random,0]);
        % Request stimulus onset
        usb.strobe;
        vbl = Screen('Flip', win, vbl + ifi/2);
        WaitSecs(stimLen);
        vbl = vbl+stimLen;
        
        Screen('DrawTexture', win,gratingTex, [],[],...
            [],[],[],[Black Black Black Black],...
            [], [],[Black,Black,Radius,centerVals(ii,1),centerVals(ii,2),checkSize,...
            random,0]);
        vbl = Screen('Flip', win, vbl + ifi/2);
        WaitSecs(waitTime);
        vbl = vbl+waitTime;
    end
end
WaitSecs(2);
usb.stopRecording;

cd('~/Documents/MATLAB/Byron/RetinoExp')
fileName = strcat('RetinoStim',Date,'_',num2str(AnimalName),'.mat');
save(fileName,'centerVals','Radius','reps','stimLen','startPause','numStimuli','w_pixels','h_pixels')
% Close window
Screen('CloseAll');

end
