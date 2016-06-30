function [] = RetinotopyCallaway(AnimalName,DistToScreen,barDegree,reps)
%RetinotopyCallaway.m
%  Display drifting horizontal and vertical bars with counter-phase
%   checkerboard patterns to map retinotopy of LFP recording electrodes or
%   single units
%  The bar will occupy 10 degrees of visual space, with an overlain checkerboard
%   pattern at 15-degree-per-side squares flashed every 167 ms, drifting so
%   that it traverses the screen in 10 seconds
%
% INPUT: AnimalName - animal's unique identifier as a number/double, e.g. 45602
%     OPTIONAL INPUTS:
%        DistToScreen - physical distance of observer from the screen, in
%           units of cm, defaults to 20 cm
%        barDegree - degrees of visual field that width of bar will
%         occupy, defaults to 10 degrees
% OUTPUT: a file named ' RetinoStim20160531_01234.mat ' where the date will be
%          adjusted for the current date and the last five digits are the
%          animal's unique ID
%
% Reference: Marshel et al. 2011 Functional specialization of seven mouse
%   visual cortical areas
%
% Created: 2016/05/31, 24 Cummington, Boston
%  Byron Price
% Updated: 2016/06/30
%  By: Byron Price

directory = pwd;
if nargin < 2
    DistToScreen = 20;
    barDegree = 10;
    reps = 20;
end
checkDegree = 15; % width or height in degrees of checkerboard squares
checkRefresh = 0.1667; % seconds to flash the checkerboard in one color
driftTime = 5;
%driftSpeed = 12; % drift speed in degrees/second

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
Width = round(((tan(barDegree*(2*pi)/360))*(DistToScreen*10))./conv_factor); % get number of pixels
                 % that barDegree degrees of visual space will occupy
                 
checkSize = round(((tan(checkDegree*(2*pi)/360))*(DistToScreen*10))./conv_factor); 

driftSpeed = [w_pixels,h_pixels]./driftTime;
%driftSpeed = ((tan(driftSpeed*(2*pi)/360))*(DistToScreen*10))./conv_factor;
                  % drift speed in pixels / second
driftSpeed = driftSpeed*ifi; % pixels / screen refresh

checkRefresh1 = round(checkRefresh/ifi);

dgshader = [directory '/RetinotopyCallaway.vert.txt'];
GratingShader = LoadGLSLProgramFromFiles({ dgshader, [directory '/RetinotopyCallaway.frag.txt'] }, 1);
gratingTex = Screen('SetOpenGLTexture', win, [], 0, GL.TEXTURE_3D,w_pixels,...
    h_pixels, 1, GratingShader);

% Define first and second ring color as RGBA vector with normalized color
% component range between 0.0 and 1.0, based on Contrast between 0 and 1
% create all textures in the same window (win), each of the appropriate
% size
Color = [0,0,0,0;1,1,1,1];

% define center positions for bar at each screen refresh
centerPos = cell(1,4);
centerPos{1} = 1:driftSpeed(1):w_pixels;
centerPos{2} = w_pixels:-driftSpeed(1):1;
centerPos{3} = 1:driftSpeed(2):h_pixels;
centerPos{4} = h_pixels:-driftSpeed(2):1;

usb.startRecording;
WaitSecs(5);

% Animation loop

% Perform initial flip to gray background and sync us to the retrace:
vbl = Screen('Flip', win);
for zz = 1:4
    centers = centerPos{zz};
    if zz == 1 || zz == 2
        vertOhorz = 1;
    else
        vertOhorz = 2;
    end
    
    for ii=1:reps
      count = 1;
      usb.strobe;
      for jj=centers
        if mod(count,checkRefresh1) <= checkRefresh1/2
            value = 1;
        else 
            value = 0;
        end
        
        Screen('DrawTexture', win,gratingTex, [],[],...
            [],[],[],[0.5 0.5 0.5 0.5],...
            [], [],[Color(1,1),Color(1,2),Color(1,3),Color(1,4),...
            Color(2,1),Color(2,2),Color(2,3),Color(2,4),Width,jj,vertOhorz,checkSize, ...
            value,0,0,0]);
            % Request stimulus onset
            vbl = Screen('Flip', win, vbl + ifi/2);
            count = count+1;
      end
      vbl = Screen('Flip',win);
      WaitSecs(1);
      vbl = vbl+1;
    end
    WaitSecs(1);
    vbl = vbl+1;
end
WaitSecs(5);
usb.stopRecording;
driftSpeed = driftSpeed/ifi; % back to pixels/second for saving purposes
stimFreq = 1/driftTime;

cd('/home/jglab/Documents/MATLAB/Byron/RetinoExp/')
fileName = strcat('RetinoStim',Date,'_',num2str(AnimalName),'.mat');
save(fileName,'driftSpeed','driftTime','stimFreq','Width','w_pixels','h_pixels','reps','checkRefresh')

% Close window
Screen('CloseAll');

end
