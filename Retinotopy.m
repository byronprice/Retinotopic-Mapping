function [] = Retinotopy(AnimalName,DistToScreen,degreeRadius)
%Retinotopy.m
%  Display a series of flashing circles to determine retinotopy of
%   LFP recording electrode.
%  Each circle will occupy a 4-degree radius of visual space
% INPUT: DistToScreen - physical distance of observer from the screen, in
%           units of cm
%        AnimalName - animal's unique identifier as a number, e.g. 45602
%        degreeRadius - degrees of visual field that radius of circle will occupy
%
% Created: 2016/05/24 at 24 Cummington, Boston
%  Byron Price
% Updated: 2016/06/30
%  By: Byron Price

directory = pwd;
if nargin < 2
    DistToScreen = 25;
    degreeRadius = 5;
    reps = 10;
    stimLen = 50/1000;
    waitTime = 1;
    startPause = 30; % 30 seconds of silence before commencing
end

Date = datetime('today','Format','yyyy-MM-dd');
Date = char(Date); Date = strrep(Date,'-','');
% Acquire a handle to OpenGL, so we can use OpenGL commands in our code:
global GL;

% Make sure this is running on OpenGL Psychtoolbox:
AssertOpenGL;

usb = usb1208FSPlusClass
WaitSecs(5);

% Choose screen with maximum id - the secondary display:
screenid = max(Screen('Screens'));

% Open a fullscreen onscreen window on that display, choose a background
% color of 128 = gray with 50% max intensity:
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
Radius = ((tan(degreeRadius*(2*pi)/360))*(DistToScreen*10))./conv_factor; % get number of pixels
     % that degreeRadius degrees of visual space will occupy
Radius = round(Radius);
centerX = Radius+1:Radius:w_pixels-Radius-1;
centerY = Radius+1:Radius:h_pixels-Radius-1;

xIndeces = randperm(length(centerX));
yIndeces = randperm(length(centerY));

numStimuli = min(length(centerX),length(centerY));
centerX = centerX(xIndeces(1:numStimuli));
centerY = centerY(yIndeces(1:numStimuli));
centerVals = zeros(length(centerX)*length(centerY),2);

count = 1;
for ii=1:length(centerX)
    for jj=1:length(centerY)
        centerVals(count,1) = centerX(ii);
        centerVals(count,2) = centerY(jj);
        count = count+1;
    end
end
permIndeces = randperm(length(centerVals));
centerVals = centerVals(permIndeces,:);
numStimuli = length(centerVals);

dgshader = [directory '/Retinotopy.vert.txt'];
GratingShader = LoadGLSLProgramFromFiles({ dgshader, [directory '/Retinotopy.frag.txt'] }, 1);
gratingTex = Screen('SetOpenGLTexture', win, [], 0, GL.TEXTURE_3D,w_pixels,...
    h_pixels, 1, GratingShader);

% Define first and second ring color as RGBA vector with normalized color
% component range between 0.0 and 1.0, based on Contrast between 0 and 1
% create all textures in the same window (win), each of the appropriate
% size
Color = [0,0,0,0;1,1,1,1];


% Perform initial flip to gray background and sync us to the retrace:
vbl = Screen('Flip', win);

usb.startRecording;
WaitSecs(startPause);

% Animation loop
for zz = 1:reps
    for ii=1:length(centerVals)
        
        % Draw the procedural texture as any other texture via 'DrawTexture'
        value = 2;
        Screen('DrawTexture', win,gratingTex, [],[],...
            [],[],[],[0 0 0 0],...
            [], [],[Color(value,1),Color(value,2),Color(value,3),Color(value,4),...
            Radius,centerVals(ii,1),centerVals(ii,2),0]);
        % Request stimulus onset
        usb.strobe;
        vbl = Screen('Flip', win, vbl + ifi/2);
        WaitSecs(stimLen);
        vbl = vbl+stimLen;
        
        Screen('DrawTexture', win,gratingTex, [],[],...
            [],[],[],[0 0 0 0],...
            [], [],[Color(1,1),Color(1,2),Color(1,3),Color(1,4),...
            Radius,centerVals(ii,1),centerVals(ii,2),0]);
        vbl = Screen('Flip', win, vbl + ifi/2);
        WaitSecs(waitTime);
        vbl = vbl+waitTime;
    end
end
WaitSecs(5);
usb.stopRecording;

fileName = strcat('RetinoStim',Date,'_',num2str(AnimalName),'.mat');
save(fileName,'centerVals','Radius','reps','stimLen','startPause','numStimuli','w_pixels','h_pixels')
% Close window
Screen('CloseAll');

end
