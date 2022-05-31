%% Thesis - Charalampos Lamprou 9114 & Ioannis Ziogas 9132 - AUTh ECE
% Thesis Project: Classification and Characterization of the Effect of Migraine 
% through Functional Connectivity Characteristics: Application to EEG 
% Recordings from a Multimorbid Clinical Sample

function [newVectors, meanValue] = remmean(vectors)

% This function removes the mean from vectors
%% Inputs: 
% vectors     -double matrix. The vectors.
%% Outputs: 
% newVectors  -double matrix. The new vectors after the substruction of the
%              mean
% meanValue   -double array. Array with the mean value of each vector
%
% This function is needed by FASTICA and FASTICAG

%-----------------------------------------------------------------------------------------------------------------
% Authors: Ioannis Ziogas & Charalampos Lamprou
% Copyright (C) 2022 Ioannis Ziogas and Charalampos Lamprou,SPBTU,ECE,AUTh
%-----------------------------------------------------------------------------------------------------------------

newVectors = zeros (size (vectors));
meanValue = mean (vectors')';
newVectors = vectors - meanValue * ones (1,size (vectors, 2));
