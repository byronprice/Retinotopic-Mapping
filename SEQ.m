function [] = SEQ(AnimalName,Day,Condition)
%SEQ.m
%  Run sequence learning
% INPUT: Obligatory-
%        AnimalName - animal's unique identifier as a number, e.g. 45602
%        Day - experimental day
%        Condition - 1 or 2 ... chosen before the experiment
%
% OUTPUT: a file with stimulus parameters named RestrictSEQStimDate_AnimalName
%           e.g. SEQStim20160708_12345.mat to be saved on the
%           CloudStation
% Created: 2017/10/02 at 24 Cummington, Boston
%  Byron Price
% Updated: 2017/10/02
%  By: Byron Price

cd('~/CloudStation/ByronExp/SEQ');
load('SEQVars.mat');

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
  
Radius = (tan((degreeRadius/2)*pi/180)*(DistToScreen*10*2))*conv_factor;

centerPosition = [w_pixels/2,h_pixels/2];

temp = (tan(((1/spatFreq)/2)*pi/180)*(DistToScreen*10*2))*conv_factor;
newSpatFreq = 1/temp;


if Day<4
    estimatedTime = ((numElements*stimTime+ISI)*reps*blocks+blocks*holdTime)/60;
    fprintf('\nEstimated time: %3.2f minutes\n',estimatedTime);
    
    % Define first and second ring color as RGBA vector with normalized color
    % component range between 0.0 and 1.0, based on Contrast between 0 and 1
    % create all textures in the same window (win), each of the appropriate
    % size
    Grey = 0.5;
    Black = 0;
    White = 1;
    
    phase = pi.*ones(numStimuli,1);
    
    stimNum = repmat([1,2,3,4],[numStimuli,1]);
    orientations = repmat([orientation(1),orientation(2),orientation(3),orientation(4)],...
        [numStimuli,1]);
    
    if Condition == 1
        interStimPause = ISI.*ones(numStimuli,1);
    elseif Condition == 2
        interStimPause = random('Uniform',ISI-0.5,ISI+0.5,[numStimuli,1]);
    else
        fprintf('\n\nCondition must be 1 or 2.\n');
        return;
    end
    
    Screen('BlendFunction',win,GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
    
    % Perform initial flip to gray background and sync us to the retrace:
    Priority(9);
    
    usb.startRecording;WaitSecs(1);usb.strobeEventWord(0);
    WaitSecs(holdTime);
    
    %Animation loop
    count = 1;
    vbl = Screen('Flip',win);
    for yy = 1:blocks
        ii=1;
        vbl = Screen('Flip',win,vbl+ifi/2);
        while ii<=reps
            
           % Draw the procedural texture as any other texture via 'DrawTexture'
            Screen('DrawTexture', win,gratingTex, [],[],...
                [],[],[],[Grey Grey Grey Grey],...
                [], [],[White,Black,...
                Radius,centerPosition(1),centerPosition(2),...
                newSpatFreq,orientations(count,1),phase(count)]);
            %Request stimulus onset
            vbl = Screen('Flip',win,vbl+ifi/2);
            %immediately strobe after stimulus onset
            usb.strobeEventWord(stimNum(count,1));
            
            Screen('DrawTexture', win,gratingTex, [],[],...
                [],[],[],[Grey Grey Grey Grey],...
                [], [],[White,Black,...
                Radius,centerPosition(1),centerPosition(2),...
                newSpatFreq,orientations(count,2),phase(count)]);
            %Request stimulus onset
            vbl = Screen('Flip',win,vbl-ifi/2+stimTime);
            usb.strobeEventWord(stimNum(count,2));
            
            Screen('DrawTexture', win,gratingTex, [],[],...
                [],[],[],[Grey Grey Grey Grey],...
                [], [],[White,Black,...
                Radius,centerPosition(1),centerPosition(2),...
                newSpatFreq,orientations(count,3),phase(count)]);
            %Request stimulus onset
            vbl = Screen('Flip',win,vbl-ifi/2+stimTime);
            usb.strobeEventWord(stimNum(count,3));
            
            Screen('DrawTexture', win,gratingTex, [],[],...
                [],[],[],[Grey Grey Grey Grey],...
                [], [],[White,Black,...
                Radius,centerPosition(1),centerPosition(2),...
                newSpatFreq,orientations(count,4),phase(count)]);
            %Request stimulus onset
            vbl = Screen('Flip',win,vbl-ifi/2+stimTime);
            usb.strobeEventWord(stimNum(count,4));
            
            vbl = Screen('Flip',win,vbl-ifi/2+stimTime);
            vbl = Screen('Flip',win,vbl-ifi/2+interStimPause(count));
            
            count = count+1;
            ii=ii+1;
            
        end
        usb.strobeEventWord(0);
        vbl = Screen('Flip',win,vbl-ifi/2+holdTime);
    end
    WaitSecs(1);
    usb.stopRecording;
    Priority(0);
    
    cd('~/CloudStation/ByronExp/SEQ');
    DayType = 'train';
    fileName = sprintf('SEQStimDay%d_%d.mat',Day,AnimalName);
    save(fileName,'spatFreq','mmPerPixel','DistToScreen','orientations',...
        'w_pixels','h_pixels','stimTime','holdTime','numElements',...
        'numStimuli','phase','stimNum','Date','DayType','ISI',...
        'interStimPause','Condition')
    % Close window
    Screen('CloseAll');

elseif Day == 4
    
    % Define first and second ring color as RGBA vector with normalized color
    % component range between 0.0 and 1.0, based on Contrast between 0 and 1
    % create all textures in the same window (win), each of the appropriate
    % size
    Grey = 0.5;
    Black = 0;
    White = 1;
    
    estimatedTime = (numConditions*(numElements*(stimTime+testStimTime)/2+ISI)*reps*blocks+numConditions*blocks*holdTime)/60;
    fprintf('\nEstimated time: %3.2f minutes\n',estimatedTime);
    
    
    orientations = repmat([orientation(1),orientation(2),orientation(3),orientation(4)],...
        [numStimuli*numConditions,1]);
    
    interStimPause = ISI.*ones(numStimuli*numConditions,1);
    
    phase = pi.*ones(numConditions*numStimuli,1);
    
    stimNum = zeros(numConditions*numStimuli,numElements);
    stimVals = [1,2,3,4;5,6,7,8;9,10,11,12;13,14,15,16];
    order = randperm(numConditions);
    
    stimTimes = stimTime*ones(numConditions*numStimuli,1);
    
    for ii=1:numConditions
        stimNum(1+(ii-1)*numStimuli:numStimuli+(ii-1)*numStimuli,:) = repmat(stimVals(order(ii),:),[numStimuli,1]);
        
        if order(ii) == 1
            interStimPause(1+(ii-1)*numStimuli:numStimuli+(ii-1)*numStimuli,:) = ISI;
        elseif order(ii) == 2
            interStimPause(1+(ii-1)*numStimuli:numStimuli+(ii-1)*numStimuli,:) = ISI;
            stimTimes(1+(ii-1)*numStimuli:numStimuli+(ii-1)*numStimuli) = testStimTime;
        elseif order(ii) == 3
            interStimPause(1+(ii-1)*numStimuli:numStimuli+(ii-1)*numStimuli,:) = ...
                random('Uniform',ISI-0.5,ISI+0.5,[numStimuli,1]);
        elseif order(ii) == 4
            interStimPause(1+(ii-1)*numStimuli:numStimuli+(ii-1)*numStimuli,:) = ...
                random('Uniform',ISI-0.5,ISI+0.5,[numStimuli,1]);
            stimTimes(1+(ii-1)*numStimuli:numStimuli+(ii-1)*numStimuli) = testStimTime;
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
        for yy = 1:blocks
            ii=1;
            vbl = Screen('Flip',win,vbl+ifi/2);
            while ii<=reps
                
                % Draw the procedural texture as any other texture via 'DrawTexture'
                Screen('DrawTexture', win,gratingTex, [],[],...
                    [],[],[],[Grey Grey Grey Grey],...
                    [], [],[White,Black,...
                    Radius,centerPosition(1),centerPosition(2),...
                    newSpatFreq,orientations(count,1),phase(count)]);
                % Request stimulus onset
                vbl = Screen('Flip',win,vbl+ifi/2);
                % immediately strobe after stimulus onset
                usb.strobeEventWord(stimNum(count,1));
                
                Screen('DrawTexture', win,gratingTex, [],[],...
                    [],[],[],[Grey Grey Grey Grey],...
                    [], [],[White,Black,...
                    Radius,centerPosition(1),centerPosition(2),...
                    newSpatFreq,orientations(count,2),phase(count)]);
                % Request stimulus onset
                vbl = Screen('Flip',win,vbl-ifi/2+stimTimes(count));
                usb.strobeEventWord(stimNum(count,2));
                
                Screen('DrawTexture', win,gratingTex, [],[],...
                    [],[],[],[Grey Grey Grey Grey],...
                    [], [],[White,Black,...
                    Radius,centerPosition(1),centerPosition(2),...
                    newSpatFreq,orientations(count,3),phase(count)]);
                % Request stimulus onset
                vbl = Screen('Flip',win,vbl-ifi/2+stimTimes(count));
                usb.strobeEventWord(stimNum(count,3));
                
                Screen('DrawTexture', win,gratingTex, [],[],...
                    [],[],[],[Grey Grey Grey Grey],...
                    [], [],[White,Black,...
                    Radius,centerPosition(1),centerPosition(2),...
                    newSpatFreq,orientations(count,4),phase(count)]);
                % Request stimulus onset
                vbl = Screen('Flip',win,vbl-ifi/2+stimTimes(count));
                usb.strobeEventWord(stimNum(count,4));
                
                vbl = Screen('Flip',win,vbl-ifi/2+stimTimes(count));
                vbl = Screen('Flip',win,vbl-ifi/2+interStimPause(count));
                
                count = count+1;
                ii=ii+1;
                
            end
            usb.strobeEventWord(0);
            vbl = Screen('Flip',win,vbl-ifi/2+holdTime);
        end
    end
    WaitSecs(1);
    usb.stopRecording;
    Priority(0);
    
    cd('~/CloudStation/ByronExp/SEQ');
    DayType = 'test';
    fileName = sprintf('SEQStimDay%d_%d.mat',Day,AnimalName);
    save(fileName,'spatFreq','mmPerPixel','DistToScreen','orientations',...
        'w_pixels','h_pixels','stimTimes','holdTime','numStimuli',...
        'phase','stimNum','Date','DayType','order','numConditions',...
        'interStimPause','Condition');
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

