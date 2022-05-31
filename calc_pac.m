%% Thesis - Charalampos Lamprou 9114 & Ioannis Ziogas 9132 - AUTh ECE
% Thesis Project: Classification and Characterization of the Effect of Migraine 
% through Functional Connectivity Characteristics: Application to EEG 
% Recordings from a Multimorbid Clinical Sample

function pac = calc_pac(param_struct)
% This function computes the phase amplitude coupling, given a Phase and an
% Amplitude signal.
% Choice over three different PAC measures is included. 
% MVL, MI and GLM measures are available. 
% NoSurrogate, Constant Thresholding, Block-Swapping and Block-Resampling
% methods are available to test statistical significance.
%
%% Inputs:
% param_struct   -struct. Contains info about pac calculation. For more 
%                 info check ConnectivityAnalysis.m
%% Outputs:
% pac            -double. The calculated pac value

%-----------------------------------------------------------------------------------------------------------------
% Authors: Ioannis Ziogas & Charalampos Lamprou
% Copyright (C) 2022 Ioannis Ziogas and Charalampos Lamprou,SPBTU,ECE,AUTh
%-----------------------------------------------------------------------------------------------------------------


Phase = param_struct.Phase;
Amp = param_struct.Amp;

if param_struct.measure == "MVL"
    z1 = (exp(1i*Phase));
    z=Amp.*(z1);% Get complex valued signal
    pac = abs((mean(z)));
elseif param_struct.measure == "MI"
    nBins = param_struct.nBins;
    [pac,~] = modulationIndex(Phase,Amp,nBins);
elseif param_struct.measure == "GLM"
    flag_AIC = param_struct.flag_AIC;
    [pac,~,~] = GLM_CFC(Phase,Amp,flag_AIC,0,0);
end

end
