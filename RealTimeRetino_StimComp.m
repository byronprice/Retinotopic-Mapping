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

numChans = 2;

tcpipClient = tcpip('128.197.59.169',30000,'NetworkRole','client');
bufferSize = 50000; % bytes, a big number (won't need this much)
set(tcpipClient,'InputBufferSize',bufferSize);
set(tcpipClient,'Timeout',5);
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

degreeRadii = [40,30,20,10,5];
degreeSpatFreq = 0.05;

% perform unit conversions
Radii = (tan(degreeRadii*pi/180)*(DistToScreen*10))*conv_factor; % get number of pixels
     % that degreeRadius degrees of visual space will occupy
temp = (tan((1/degreeSpatFreq)*pi/180)*(DistToScreen*10))*conv_factor;
spatFreq = 1/temp;

numStimuli = 4; % must be a perfect square
centerVals = zeros(numStimuli,2);
xaxis = w_pixels/(numStimuli-1);
yaxis = h_pixels/(numStimuli-1);

count = 1;
for ii=1:sqrt(numStimuli)
    for jj=1:sqrt(numStimuli)
        centerVals(count,1) = xaxis*ii;
        centerVals(count,2) = yaxis*jj;
        count = count+1;
    end
end

Data = zeros(numStimuli,1);
WaitSecs(15);

usb.strobeEventWord(startEXP);
Pr = fread(tcpipClient,numChans,'double');

Priority(9);
bigData = [];
for ii=1:numChans
    for jj=1:10
        usb.strobeEventWord(jj);
        dataSize = fread(tcpipClient,3,'double');
        
        if isempty(dataSize) == 0
            temp = fread(tcpipClient,dataSize(1)*dataSize(2),'double');
            data = reshape(temp,[dataSize(1),dataSize(2)]);
            bigData = [bigData;data];
        end
    end
    usb.strobeEventWord(endRUN);
end

Priority(0);
Screen('CloseAll');
cd('~/CloudStation/ByronExp/Retino');
Date = datetime('today','Format','yyyy-MM-dd');
Date = char(Date); Date = strrep(Date,'-','');Date=str2double(Date);

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