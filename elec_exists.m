%% Thesis - Charalampos Lamprou 9114 & Ioannis Ziogas 9132 - AUTh ECE
% Thesis Project: Classification and Characterization of the Effect of Migraine 
% through Functional Connectivity Characteristics: Application to EEG 
% Recordings from a Multimorbid Clinical Sample

function [flag,pos] = elec_exists(elec_mat, elec)
% This function checks if a specific electrode exists in a string array
%% Inputs:
% elec_mat    -string array. Contains the names of all possible electrodes
% elec        -string. A specific electrode
%% Outputs:
% flag        -binary flag. 1 if the elec is in the elec_mat and 0
%              otherwise
% pos         -double. The position of elec in elec_mat
%
%-----------------------------------------------------------------------------------------------------------------
% Authors: Ioannis Ziogas & Charalampos Lamprou
% Copyright (C) 2022 Ioannis Ziogas and Charalampos Lamprou,SPBTU,ECE,AUTh
%-----------------------------------------------------------------------------------------------------------------

flag = 0;
counter = 1;
pos = NaN;
for i = 1:length(elec_mat)
    if elec_mat{i} == elec
        flag = 1;
        pos = counter;
    end
    counter = counter + 1;
end