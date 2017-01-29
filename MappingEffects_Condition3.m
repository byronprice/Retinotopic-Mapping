function [] = MappingEffects_Condition3(AnimalName,ConditionNumber)
%MappingEffects_Condition3.m
%  Third condition for MappingEffects experiment
%   Fake mapping stimulus, followed by SRP protocol
cd('~/CloudStation/ByronExp/MappingEffects');
load('RetinotopyVars.mat');
load('SRPvars.mat');

directory = '~/Documents/MATLAB/Byron/Retinotopic-Mapping/';

holdTime = 30; % 30 second pauses between blocks

Date = datetime('today','Format','yyyy-MM-dd');
Date = char(Date); Date = strrep(Date,'-','');Date=str2double(Date);
% Acquire a handle to OpenGL, so we can use OpenGL commands in our code:
global GL;

% Make sure this is running on OpenGL Psychtoolbox:
AssertOpenGL;

% usb = ttlInterfaceClass.getTTLInterface;
usb = usb1208FSPlusClass;
display(usb);

WaitSecs(10);
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

dgshader = [directory '/SRP.vert.txt'];
GratingShader = LoadGLSLProgramFromFiles({ dgshader, [directory '/SRP.frag.txt'] }, 1);
gratingTex = Screen('SetOpenGLTexture', win, [], 0, GL.TEXTURE_3D,w_pixels,...
    h_pixels, 1, GratingShader);

% screen size in millimeters and a conversion factor to get from mm to pixels
[w_mm,h_mm] = Screen('DisplaySize',screenid);
conv_factor = (w_mm/w_pixels+h_mm/h_pixels)/2;
mmPerPixel = conv_factor;
conv_factor = 1/conv_factor;

% perform unit conversions
temp = (tan((1/spatFreq)*pi/180)*(DistToScreen*10))*conv_factor;
spatFreq = 1/temp;

% Define first and second ring color as RGBA vector with normalized color
% component range between 0.0 and 1.0, based on Contrast between 0 and 1
% create all textures in the same window (win), each of the appropriate
% size
Grey = 0.5;
Black = 0;
White = 1;

% Perform initial flip to gray background and sync us to the retrace:
Priority(9);

usb.startRecording;WaitSecs(1);usb.strobeEventWord(0);
WaitSecs(holdTime);

WaitSecs(30*60);

vbl = Screen('Flip',win);
for yy = 1:srp_blocks
    phase = rand*pi;
    for zz = 1:srp_reps/srp_blocks
        % Draw the procedural texture as any other texture via 'DrawTexture'
        Screen('DrawTexture', win,gratingTex, [],[],...
            [],[],[],[Grey Grey Grey Grey],...
            [], [],[White,Black,...
            spatFreq,srp_orientation,phase,0,0,0]);
        % Request stimulus onset
        vbl = Screen('Flip', win,vbl-ifi/2+1/srp_Hz);usb.strobeEventWord(100);
        phase = phase+pi;
    end
    vbl = Screen('Flip',win,vbl+ifi/2);usb.strobeEventWord(0);
    vbl = Screen('Flip',win,vbl-ifi/2+holdTime);
end
WaitSecs(2);
usb.stopRecording;
Priority(0);

fileName = sprintf('MappingEffectsStim%d_%d.mat',Date,AnimalName);
save(fileName,'ConditionNumber','srp_blocks','srp_orientation','srp_reps','spatFreq',...
    'srp_Hz');
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
