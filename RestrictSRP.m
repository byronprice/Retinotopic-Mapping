function [] = RestrictSRP(AnimalName,Day)
%RestrictSRP.m
%  Run SRP within a restricted region surrounding the LFP retinotopic response
%    region.
% INPUT: Obligatory-
%        AnimalName - animal's unique identifier as a number, e.g. 45602
%        Day - experimental day
%
% OUTPUT: a file with stimulus parameters named RestrictSRPStimDate_AnimalName
%           e.g. RestrictSRPStim20160708_12345.mat to be saved on the
%           CloudStation
% Created: 2017/08/08 at 24 Cummington, Boston
%  Byron Price
% Updated: 2017/08/08
%  By: Byron Price

cd('~/CloudStation/ByronExp/RestrictSRP');
load('RestrictSRPVars.mat');

currentdirectory = '~/Documents/MATLAB/Byron/Retinotopic-Mapping';
cd(currentdirectory);

reps = numStimuli/blocks;

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

if Day==1
    [centerPositions,targetChan] = GetRetinoMap(AnimalName);
    targetChan = 1;
else
   cd('~/CloudStation/ByronExp/RestrictSRP');
   fileName = sprintf('RestrictSRPStimDay1_%d.mat',AnimalName);
   load(fileName,'centerPositions','targetChan');
   cd(currentdirectory);
end


if Day<5
    estimatedTime = (stimTime*reps*blocks+blocks*holdTime)/60;
    fprintf('\nEstimated time: %3.2f minutes\n',estimatedTime);
    
    % Define first and second ring color as RGBA vector with normalized color
    % component range between 0.0 and 1.0, based on Contrast between 0 and 1
    % create all textures in the same window (win), each of the appropriate
    % size
    Grey = 0.5;
    Black = 0;
    White = 1;
    
    phase = pi.*ones(numStimuli,1);
    phase(1:2:end) = 0;
    
    stimNum = 2.*ones(numStimuli,1);
    stimNum(1:2:end) = 1;
    
    Screen('BlendFunction',win,GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
    
    % Perform initial flip to gray background and sync us to the retrace:
    Priority(9);
    
    usb.startRecording;WaitSecs(1);usb.strobeEventWord(0);
    WaitSecs(holdTime);
    
    % Animation loop
    count = 1;
    vbl = Screen('Flip',win);
    for yy = 1:blocks
        vbl = Screen('Flip',win,vbl+ifi/2);
        ii=1;
        while ii<=reps
            
            % Draw the procedural texture as any other texture via 'DrawTexture'
            Screen('DrawTexture', win,gratingTex, [],[],...
                [],[],[],[Grey Grey Grey Grey],...
                [], [],[White,Black,...
                Radius,centerPositions(targetChan,1),centerPositions(targetChan,2),newSpatFreq,orientation,phase(count)]);
            % Request stimulus onset
            vbl = Screen('Flip',win,vbl-ifi/2+stimTime);
            usb.strobeEventWord(stimNum(count));
            count = count+1;
            ii=ii+1;
            
        end
        vbl = Screen('Flip',win,vbl-ifi/2+stimTime);
        usb.strobeEventWord(0);
        vbl = Screen('Flip',win,vbl-ifi/2+holdTime);
    end
    WaitSecs(1);
    usb.stopRecording;
    Priority(0);
    
    cd('~/CloudStation/ByronExp/RestrictSRP');
    DayType = 'train';
    fileName = sprintf('RestrictSRPStimDay%d_%d.mat',Day,AnimalName);
    save(fileName,'centerPositions','targetChan','Radius','degreeRadius','spatFreq',...
        'mmPerPixel','DistToScreen','orientation','w_pixels','h_pixels','stimTime','holdTime',...
        'numStimuli','phase','stimNum','Date','DayType')
    % Close window
    Screen('CloseAll');

elseif Day == 5
    
    % Define first and second ring color as RGBA vector with normalized color
    % component range between 0.0 and 1.0, based on Contrast between 0 and 1
    % create all textures in the same window (win), each of the appropriate
    % size
    Grey = 0.5;
    Black = 0;
    White = 1;
    
    numConditions = 4;
    
    estimatedTime = (numConditions*stimTime*reps*blocks+numConditions*blocks*holdTime)/60;
    fprintf('\nEstimated time: %3.2f minutes\n',estimatedTime);
    
    
    phase = pi.*ones(numConditions*numStimuli,1);
    phase(1:2:end) = 0;
    
    stimNum = zeros(numConditions*numStimuli,1);
    stimVals = [1,2;3,4;5,6;7,8];
    order = randperm(numConditions);
    channel = zeros(numConditions,1);
    
    orientations = orientation.*ones(4*numStimuli,1);
    for ii=1:numConditions
        stimNum(1+(ii-1)*numStimuli:2:numStimuli+(ii-1)*numStimuli) = stimVals(order(ii),1);
        stimNum(2+(ii-1)*numStimuli:2:numStimuli+(ii-1)*numStimuli) = stimVals(order(ii),2);
        
        if order(ii) == 1
            channel(ii) = targetChan;
        elseif order(ii) == 2
            channel(ii) = targetChan;
            orientations(1+(ii-1)*numStimuli:numStimuli+(ii-1)*numStimuli) = orientation+3*pi/4;
        elseif order(ii) == 3
            channel(ii) = -targetChan+3;
        elseif order(ii) == 4
            channel(ii) = -targetChan+3;
            orientations(1+(ii-1)*numStimuli:numStimuli+(ii-1)*numStimuli) = orientation+pi/4;
        end
    end
    
    
    Screen('BlendFunction',win,GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
    
    % Perform initial flip to gray background and sync us to the retrace:
    Priority(9);
    
    usb.startRecording;WaitSecs(1);usb.strobeEventWord(0);
    WaitSecs(holdTime);
    
    % Animation loop
    count = 1;
    vbl = Screen('Flip',win);
    for zz=1:numConditions
        vbl = Screen('Flip',win,vbl+ifi/2);
        for yy = 1:blocks
            ii=1;
            while ii<=reps
                
                % Draw the procedural texture as any other texture via 'DrawTexture'
                Screen('DrawTexture', win,gratingTex, [],[],...
                    [],[],[],[Grey Grey Grey Grey],...
                    [], [],[White,Black,...
                    Radius,centerPositions(channel(zz),1),centerPositions(channel(zz),2),...
                    newSpatFreq,orientations(count),phase(count)]);
                % Request stimulus onset
                vbl = Screen('Flip',win,vbl-ifi/2+stimTime);
                usb.strobeEventWord(stimNum(count));
                count = count+1;
                ii=ii+1;
                
            end
            vbl = Screen('Flip',win,vbl-ifi/2+stimTime);
            usb.strobeEventWord(0);
            vbl = Screen('Flip',win,vbl-ifi/2+holdTime);
        end
    end
    WaitSecs(1);
    usb.stopRecording;
    Priority(0);
    
    cd('~/CloudStation/ByronExp/RestrictSRP');
    DayType = 'test';
    fileName = sprintf('RestrictSRPStimDay%d_%d.mat',Day,AnimalName);
    save(fileName,'centerPositions','targetChan','Radius','degreeRadius','spatFreq',...
        'mmPerPixel','DistToScreen','orientations','w_pixels','h_pixels','stimTime','holdTime',...
        'numStimuli','phase','stimNum','Date','DayType','order','channel','numConditions')
    % Close window
    Screen('CloseAll');
    
end

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

function [centerPositions,targetChan] = GetRetinoMap(AnimalName)
cd ~/CloudStation/ByronExp/Retino/
fileName = strcat('RetinoMapBayes*',num2str(AnimalName),'.mat');
files = dir(fileName);

load(files(end).name);

[numChans,~,numSamples] = size(posteriorSample);

x = 1:w_pixels;y=1:h_pixels;
[X,Y] = meshgrid(x,y); 

% we want final RF to be ~1000 by 1000 screen pixels, about 50 degrees
%  of visual arc on a side

centerPositions = zeros(numChans,2);
bestChan = zeros(numChans,1);
for ii=1:numChans
    finalIm = zeros(length(y),length(x));
    samples = squeeze(posteriorSample(ii,:,:));
    N = 1000;
    for ll=1:N
        index = random('Discrete Uniform',numSamples);
        parameterVec = samples(:,index);
        b = [parameterVec(1),parameterVec(4),parameterVec(5),parameterVec(6)];
        distX = X-parameterVec(2);distY = Y-parameterVec(3);
        finalIm = finalIm+b(1)*exp(-(distX.^2)./(2*b(2)*b(2))-...
                     (distY.^2)./(2*b(3)*b(3)))+b(4);
%         for jj=1:length(x)
%             for kk=1:length(y)
%                 distX = x(jj)-parameterVec(2);
%                 distY = y(kk)-parameterVec(3);
%                 
%                 finalIm(jj,kk) = finalIm(jj,kk)+b(1)*exp(-(distX.^2)./(2*b(2)*b(2))-...
%                     (distY.^2)./(2*b(3)*b(3)))+b(4);
%             end
%         end
    end
    finalIm = finalIm./N;
    [~,maxInd] = max(finalIm(:));
    [row,col] = ind2sub(size(finalIm),maxInd);
    centerPositions(ii,1) = col;
    centerPositions(ii,2) = row;
    
    [~,minInd] = min(finalIm(:));
    
    bestChan(ii) = finalIm(maxInd)-finalIm(minInd);
end

[~,targetChan] = max(bestChan);

cd ~/CloudStation/ByronExp/RestrictSRP/
end
