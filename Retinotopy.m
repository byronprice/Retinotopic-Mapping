function [] = Retinotopy(AnimalName,DistToScreen,degreeRadius)
%Retinotopy.m
%  Display a series of flashing circles to determine retinotopy of
%   LFP recording electrode.
%  Each circle will occupy a 1.5 degree radius of visual space
% INPUT: DistToScreen - physical distance of observer from the screen, in
%           units of cm
%        AnimalName - animal's unique identifier as a number, e.g. 45602
%        degreeRadius - degrees of visual that radius of circle will occupy
%
% Created: 2016/05/24 at 24 Cummington, Boston
%  Byron Price
% Updated: 2016/05/25
%  By: Byron Price

directory = pwd;
if nargin < 2
    DistToScreen = 20;
    degreeRadius = 2;
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
[win,~] = Screen('OpenWindow', screenid,128);

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
     % that 1 degree of visual space will occupy
Radius = round(Radius);
CenterX = w_pixels/2-400:Radius*2:w_pixels/2+400;
CenterY = h_pixels/2-400:Radius*2:h_pixels/2+400;

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
ts = vbl;

reps = 5;
stimulusLocs = zeros(length(CenterX)*length(CenterX)*reps,2);

usb.startRecording;
WaitSecs(5);

% Animation loop
count = 1;
for zz = 1:reps
    for ii=1:length(CenterX)
      for jj=1:length(CenterY)
        stimulusLocs(count,:) = [CenterX(ii),CenterY(jj)];
        count = count+1;
        % Draw the procedural texture as any other texture via 'DrawTexture'
        usb.strobe;
        for kk=1:8
            value = mod(kk,2)+1;
            Screen('DrawTexture', win,gratingTex, [],[],...
                [],[],[],[0.5 0.5 0.5 0.5],...
                [], [],[Color(value,1),Color(value,2),Color(value,3),Color(value,4),...
                Radius,CenterX(ii),CenterY(jj),0]);
                % Request stimulus onset
                vbl = Screen('Flip', win, vbl + 5*ifi/2);

                %usb.triggerON(1,7);
                %usb.triggerOFF(1,7);
        end
        usb.strobe;

         Screen('DrawTexture', win,gratingTex, [],[],...
             [],[],[],[0.5 0.5 0.5 0.5],...
             [], [],[0.5,0.5,0.5,0.5,...
             Radius,CenterX(ii),CenterY(jj),0]);
        vbl = Screen('Flip', win, vbl + 5*ifi/2);
        WaitSecs(1.0);
        vbl = vbl+1.0;
      end
    end
end
WaitSecs(5);
usb.stopRecording;

fileName = strcat('RetinoStim',Date,'_',num2str(AnimalName),'.mat');
save(fileName,'stimulusLocs','Radius')
% Close window
Screen('CloseAll');

end
