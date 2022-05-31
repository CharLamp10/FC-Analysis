%% Thesis - Charalampos Lamprou 9114 & Ioannis Ziogas 9132 - AUTh ECE
% Thesis Project: Classification and Characterization of the Effect of Migraine 
% through Functional Connectivity Characteristics: Application to EEG 
% Recordings from a Multimorbid Clinical Sample

function bands = create_bands_pre(lower_limit,upper_limit,bandwidth,overlap)
% Function to create the desired non conventional bands
% with specified bandwidth and overlap.
%% Inputs:
% lower_limit     -double. The lower frequency of interest in Hz
% upper_limit     -double. The higher frequency of interest in Hz
% bandwidth       -double. The bandwidth of each band in Hz
% overlap         -double. The overlap between bands
%% Outputs:
% bands           -double array. The limits of the non-conventional bands
% 
%-----------------------------------------------------------------------------------------------------------------
% Authors: Ioannis Ziogas & Charalampos Lamprou
% Copyright (C) 2022 Ioannis Ziogas and Charalampos Lamprou,SPBTU,ECE,AUTh
%-----------------------------------------------------------------------------------------------------------------


bands(1) = lower_limit;
bands(2) = lower_limit + bandwidth;

counter = 3;
while bands(end) < upper_limit 
    bands(counter) = bands(end) - bandwidth*(overlap/100);
    counter = counter + 1;
    bands(counter) = bands(end) + bandwidth;
    counter = counter + 1;
end

end