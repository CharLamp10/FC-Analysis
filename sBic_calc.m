%% Thesis - Charalampos Lamprou 9114 & Ioannis Ziogas 9132 - AUTh ECE
% Thesis Project: Classification and Characterization of the Effect of Migraine 
% through Functional Connectivity Characteristics: Application to EEG 
% Recordings from a Multimorbid Clinical Sample
function totalBic = sBic_calc(bic,upLim1,upLim2,lowLim1,lowLim2)
% This function calculates the total Bicoherence/PAC values in a given
% frequency block.
%% Inputs:
% bic      -double matrix. A matrix with the feature values. For each
%           frequency pair there is a feature value.
% upLim1   -double. upper frequency for the x axis
% upLim2   -double. upper frequency for the y axis
% lowLim1  -double. lower frequency for the x axis
% lowLim2  -double. lower frequency for the y axis
%% Outputs:
% totalBic -double. The total feature value in the given block.

%-----------------------------------------------------------------------------------------------------------------
% Authors: Ioannis Ziogas & Charalampos Lamprou
% Copyright (C) 2022 Ioannis Ziogas and Charalampos Lamprou,SPBTU,ECE,AUTh
%-----------------------------------------------------------------------------------------------------------------


lowLim1 = ceil(lowLim1); lowLim2 = ceil(lowLim2); 
upLim1 = ceil(upLim1); upLim2 = ceil(upLim2);
bic = bic(lowLim1:upLim1,lowLim2:upLim2);
totalBic = sum(sum(bic));

end