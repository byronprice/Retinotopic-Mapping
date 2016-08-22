function [] = SlidingBars(Length,Width,Orientation,Foreground,Background,DistToScreen)
%SlidingBars.m
%  Display drifting bars with the option to add a checkered pattern
%
% INPUT: Length - bar length (degress of visual space)
%        Width - bar width (degrees of visual space)
%        Orientation - bar orientation (degrees)
%        Foreground - color of the bar (0 to 1, 0 is black, 1 is white)
%        Background - color of the background (as above)
%        Checkered - 'yes' or 'no'
% OUTPUT: none
%
%
% Created: 2016/08/19, 24 Cummington, Boston
%  Byron Price
% Updated: 2016/08/19
%  By: Byron Price

directory = '~/Documents/Current-Projects/Retinotopic-Mapping';

% Acquire a handle to OpenGL, so we can use OpenGL commands in our code:
global GL;

if nargin < 6
    DistToScreen = 25;
end

gama = 2.1806;

% if strcmp(Checkered,'yes') == 1
%     Checkered = 1;
%     checkDegree = Width*2;
%     checkRefresh = 0.5;
% end
% Make sure this is running on OpenGL Psychtoolbox:
AssertOpenGL;

% Choose screen with maximum id - the secondary display:
screenid = max(Screen('Screens'));

% % Open a fullscreen onscreen window on that display, choose a background
% % color of 127 = gray with 50% max intensity; 0 = black;255 = white
background = round(Background*255);
[win,~] = Screen('OpenWindow', screenid,background);

gammaTable = makeGrayscaleGammaTable(gama,0,255);
Screen('LoadNormalizedGammaTable',win,gammaTable);

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
Width = round(((tan(Width*(2*pi)/360))*(DistToScreen*10))./conv_factor); % get number of pixels
                 % that Width degrees of visual space will occupy
                 
Length = round(((tan(Length*(2*pi)/360))*(DistToScreen*10))./conv_factor);
                 
%checkSize = round(((tan(checkDegree*(2*pi)/360))*(DistToScreen*10))./conv_factor); 

dgshader = [directory '/SlidingBars.vert.txt'];
GratingShader = LoadGLSLProgramFromFiles({ dgshader, [directory '/SlidingBars.frag.txt'] }, 1);
gratingTex = Screen('SetOpenGLTexture', win, [], 0, GL.TEXTURE_3D,w_pixels,...
    h_pixels, 1, GratingShader);

% Define first and second ring color as RGBA vector with normalized color
% component range between 0.0 and 1.0, based on Contrast between 0 and 1
% create all textures in the same window (win), each of the appropriate
% size

% define center positions for bar at each screen refresh
centerPos = [w_pixels/2 h_pixels/2];

Screen('BlendFunction',win,GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);

floatSpeed = 5;
target = [rand*w_pixels rand*h_pixels];
runTime = 60; % seconds
randVals = rand([round((1/(ifi/1000)*runTime)),2]);
randVals = randVals*[w_pixels,0;0,h_pixels];
% Animation loop
count = 1;
vbl = Screen('Flip', win);
while ~KbCheck
    Screen('DrawTexture', win,gratingTex, [],[],...
        [],[],[],[0.5 0.5 0.5 0.5],...
        [], [],[Foreground,Background,Width,Length,centerPos(1),centerPos(2),...
        Orientation,0]);
    % Request stimulus onset
    vbl = Screen('Flip', win,vbl+ifi/2);
    vec = target-centerPos;
    mag = norm(vec);
    centerPos = centerPos+(vec./mag)*floatSpeed;
    
    if mag <= floatSpeed/2
        target(1) = randVals(count,1);
        target(2) = randVals(count,2);
    end
    Orientation = Orientation+1;
    count = count+1;
end
WaitSecs(2);

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

