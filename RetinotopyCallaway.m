function [] = RetinotopyCallaway(AnimalName,holdTime)
%RetinotopyCallaway.m
%  Display drifting horizontal and vertical bars with counter-phase
%   checkerboard patterns to map retinotopy of LFP recording electrodes or
%   single units
%  The bar will occupy ~10 degrees of visual space, with an overlain checkerboard
%   pattern at 15-degree-per-side squares flashed every 200 ms, drifting so
%   that it traverses the screen in 10 seconds
%
% INPUT: AnimalName - animal's unique identifier as a number/double, e.g. 45602
%     OPTIONAL INPUTS:
%        holdTime - amount of time to wait at the beginning and between
%        changes in direction
%      other changeable variables in the file RetinotopyCallawayVars.mat
% OUTPUT: a file named ' RetinoStim20160531_01234.mat ' where the date will be
%          adjusted for the current date and the last five digits are the
%          animal's unique ID
%
% Reference: Marshel et al. 2011 Functional specialization of seven mouse
%   visual cortical areas
%
% Created: 2016/05/31, 24 Cummington, Boston
%  Byron Price
% Updated: 2016/08/17
%  By: Byron Price

cd('~/CloudStation/ByronExp/Retino');
load('RetinotopyCallawayVars.mat');

directory = '~/Documents/MATLAB/Byron/Retinotopic-Mapping';
%directory = '~/Documents/Current-Projects/Retinotopic-Mapping';

if nargin < 2
   holdTime = 30;
end

Date = datetime('today','Format','yyyy-MM-dd');
Date = char(Date); Date = strrep(Date,'-','');Date=str2double(Date);
% Acquire a handle to OpenGL, so we can use OpenGL commands in our code:
global GL;

% Make sure this is running on OpenGL Psychtoolbox:
AssertOpenGL;

usb = usb1208FSPlusClass;
WaitSecs(10);

% Choose screen with maximum id - the secondary display:
screenid = max(Screen('Screens'));

% % Open a fullscreen onscreen window on that display, choose a background
% % color of 127 = gray with 50% max intensity; 0 = black;255 = white
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

% screen size in millimeters and a conversion factor to get from mm to pixels
[w_mm,h_mm] = Screen('DisplaySize',screenid);
conv_factor = (w_mm/w_pixels+h_mm/h_pixels)/2;

% perform unit conversions
Width = round(((tan(barDegree*(2*pi)/360))*(DistToScreen*10))./conv_factor); % get number of pixels
                 % that barDegree degrees of visual space will occupy
                 
checkSize = round(((tan(checkDegree*(2*pi)/360))*(DistToScreen*10))./conv_factor); 

driftSpeed = [w_pixels,h_pixels]./driftTime;
%driftSpeed = ((tan(driftSpeed*(2*pi)/360))*(DistToScreen*10))./conv_factor;
                  % drift speed in pixels / second
driftSpeed = driftSpeed.*ifi; % pixels / screen refresh

checkRefresh1 = round(checkRefresh/ifi);

dgshader = [directory '/RetinotopyCallaway.vert.txt'];
GratingShader = LoadGLSLProgramFromFiles({ dgshader, [directory '/RetinotopyCallaway.frag.txt'] }, 1);
gratingTex = Screen('SetOpenGLTexture', win, [], 0, GL.TEXTURE_3D,w_pixels,...
    h_pixels, 1, GratingShader);

% Define first and second ring color as RGBA vector with normalized color
% component range between 0.0 and 1.0, based on Contrast between 0 and 1
% create all textures in the same window (win), each of the appropriate
% size
Color = [0,1];

% define center positions for bar at each screen refresh
numDirs = 4;
DirNames = {'Right','Left','Up','Down'};
centerPos = cell(numDirs,1);
centerPos{1} = 1:driftSpeed(1):w_pixels;
centerPos{2} = w_pixels:-driftSpeed(1):1;
centerPos{3} = 1:driftSpeed(2):h_pixels;
centerPos{4} = h_pixels:-driftSpeed(2):1;

Flashes = cell(numDirs,1);
for ii=1:numDirs
    centers = centerPos{ii};
    values = mod(0:length(centers)-1,checkRefresh1) < checkRefresh1/2;
    temp = diff(values);Flashes{ii} = [0,temp];
end

eventNums = cell(numDirs,1);

estimatedTime = ((driftTime+1)*reps+holdTime)*numDirs/60;
display(sprintf('\nEstimated time: %3.2f minutes',estimatedTime));

Screen('BlendFunction',win,GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);

usb.startRecording;WaitSecs(1);usb.strobeEventWord(0);
WaitSecs(holdTime);

% Animation loop
vbl = Screen('Flip', win);
for zz = 1:numDirs
    centers = centerPos{zz};
    if zz == 1 || zz == 2
        vertOhorz = 1;
    else
        vertOhorz = 2;
    end
    for ii=1:reps
        vbl = Screen('Flip',win,vbl+ifi/2);
      for jj=1:length(centers)
        Screen('DrawTexture', win,gratingTex, [],[],...
            [],[],[],[0.5 0.5 0.5 0.5],...
            [], [],[Color(1),Color(2),Width,centers(jj),vertOhorz,checkSize, ...
            values(jj),0]);
            % Request stimulus onset
            if abs(Flashes(jj)) > 0 
                usb.strobeEventWord(zz);
            end
            vbl = Screen('Flip', win,vbl+ifi/2);
      end
      vbl = Screen('Flip', win,vbl+ifi/2);
      vbl = Screen('Flip',win,vbl-ifi/2+1);
    end
    if zz ~= numDirs
        usb.strobeEventWord(0);
        vbl = Screen('Flip',win,vbl-ifi/2+holdTime);
    end
end
WaitSecs(5);
usb.stopRecording;
driftSpeed = driftSpeed/ifi; % back to pixels/second for saving purposes
stimFreq = 1/driftTime;

stimParams = RetinoCallStimObj;
stimParams.driftSpeed = driftSpeed;
stimParams.driftTime = driftTime;
stimParams.stimFreq = stimFreq;
stimParams.Width = Width;
stimParams.w_pixels = w_pixels;
stimParams.h_pixels = h_pixels;
stimParams.reps = reps;
stimParams.checkRefresh = checkRefresh;
stimParams.holdTime = holdTime;
stimParams.centerPos = centerPos;
stimParams.Flashes = Flashes;
stimParams.numDirs = numDirs;
stimParams.DistToScreen = DistToScreen;
stimParams.DirNames = DirNames;

cd('~/CloudStation/ByronExp/Retino/')
fileName = sprintf('RetinoCallStim%d_%d.mat',Date,AnimalName);
save(fileName,'stimParams')

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
