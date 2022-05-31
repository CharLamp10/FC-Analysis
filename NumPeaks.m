%% Thesis - Charalampos Lamprou 9114 & Ioannis Ziogas 9132 - AUTh ECE
% Thesis Project: Classification and Characterization of the Effect of Migraine 
% through Functional Connectivity Characteristics: Application to EEG 
% Recordings from a Multimorbid Clinical Sample

function noPeaks = NumPeaks(Peaks,upLim1,upLim2,lowLim1,lowLim2)
% This function calculates the number of peaks in a given frequency block
%% Inputs:
% Peaks     -double array. Contains the coordinates of each peak
% upLim1    -double. upper frequency for the x axis
% upLim2    -double. upper frequency for the y axis
% lowLim1   -double. lower frequency for the x axis
% lowLim2   -double. lower frequency for the y axis
%% Outputs:
% noPeaks   -double. The number of peaks

%-----------------------------------------------------------------------------------------------------------------
% Authors: Ioannis Ziogas & Charalampos Lamprou
% Copyright (C) 2022 Ioannis Ziogas and Charalampos Lamprou,SPBTU,ECE,AUTh
%-----------------------------------------------------------------------------------------------------------------

noPeaks = 0;
for i = 1:length(Peaks(1,:))
    if Peaks(1,i) > lowLim1 && Peaks(1,i) < upLim1...
            && Peaks(2,i) > lowLim2 && Peaks(2,i) < upLim2
        noPeaks = noPeaks + 1;
    end
end

end