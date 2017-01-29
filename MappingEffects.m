function [] = MappingEffects(AnimalName,ConditionNumber)
%MappingEffects.m
%  Experiment that Jeff Gavornik and I devised to measure the effects of presenting
%   the mapping stimulus (does it potentiate? ... does it affect the
%   potentiation to the SRP protocol?)
%  For the mapping stimulus, a series of circles are displayed to the mouse
%   Each circle will occupy an ~ 5-degree radius of visual space and be a
%   black/white 2D sinusoid masked by a Gaussian kernel
% INPUT: Obligatory-
%        AnimalName - animal's unique identifier as a number, e.g. 45602
%        ConditionNumber - 1, 2, or 3 for the experimental condition
%           1) just mapping followed by fake SRP (grey screen)
%           2) mapping followed by SRP
%           3) fake mapping followed by SRP
%        see file RetinotopyVars.mat for other changeable presets
%
% OUTPUT: a file with stimulus parameters named MappingEffectsStimDate_AnimalName
%           e.g. MappingEffectsStim20160708_12345.mat to be saved in the
%           MappingEffects
%           folder in the CloudStation
% Created: 2017/01/28 at 24 Cummington, Boston
%  Byron Price
% Updated: 2017/01/28
%  By: Byron Price

if ConditionNumber == 1
    MappingEffects_Condition1(AnimalName,ConditionNumber);
elseif ConditionNumber == 2
    MappingEffects_Condition2(AnimalName,ConditionNumber);
elseif ConditionNumber == 3
    MappingEffects_Condition3(AnimalName,ConditionNumber);
else
    display('ConditionNumber must be either 1,2, or 3');
end

end
