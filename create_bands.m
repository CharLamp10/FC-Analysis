%% Thesis - Charalampos Lamprou 9114 & Ioannis Ziogas 9132 - AUTh ECE
% Thesis Project: Classification and Characterization of the Effect of Migraine 
% through Functional Connectivity Characteristics: Application to EEG 
% Recordings from a Multimorbid Clinical Sample

function bands = create_bands(lower_limit,upper_limit,bandwidth,overlap,untouched)
% Function to create the desired non conventional bands
% with specified bandwidth and overlap.
%% Inputs:
% lower_limit     -double. The lower frequency of interest in Hz
% upper_limit     -double. The higher frequency of interest in Hz
% bandwidth       -double. The bandwidth of each band in Hz
% overlap         -double. The overlap between bands
% untouched       -double array. Limits in which different band parameters
%                  will be applied
%% Outputs:
% bands           -double array. The limits of the non-conventional bands
% 
%-----------------------------------------------------------------------------------------------------------------
% Authors: Ioannis Ziogas & Charalampos Lamprou
% Copyright (C) 2022 Ioannis Ziogas and Charalampos Lamprou,SPBTU,ECE,AUTh
%-----------------------------------------------------------------------------------------------------------------

if length(bandwidth) == 1
    if length(untouched) >= 1
        error('utouched should be empty')
    end
    bands = create_bands_pre(lower_limit,upper_limit,bandwidth,overlap);
else
    if length(bandwidth) - 1 ~= length(untouched)
        error('The length of the untouched must be equal to length of bandwidth - 1\n')
    end
      
    untouched(end + 1) = upper_limit;
    bands = create_bands_pre(lower_limit,untouched(1),bandwidth(1),overlap);
    for i = 2:length(bandwidth)
        temp = bands(end) - bandwidth(i)*overlap/100;
        bands = [bands create_bands_pre(temp,untouched(i),bandwidth(i),overlap)];
    end

end

if bands(end) > upper_limit && (bands(end-2) + round(0.5*bandwidth) < upper_limit)
    bands(end) = upper_limit;
elseif bands(end) == upper_limit
    
else
    bands = bands(1:end-2);
end
end