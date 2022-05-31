%% Thesis - Charalampos Lamprou 9114 & Ioannis Ziogas 9132 - AUTh ECE
% Thesis Project: Classification and Characterization of the Effect of Migraine 
% through Functional Connectivity Characteristics: Application to EEG 
% Recordings from a Multimorbid Clinical Sample

function Ent = bspecEntropy(mat,upLim1,upLim2,lowLim1,lowLim2,q)
% Calculates the Entropy in a specific frequency block
%% Inputs: 
% mat      -double matrix. The given feature matrix
% upLim1   -double. upper frequency for the x axis
% upLim2   -double. upper frequency for the y axis
% lowLim1  -double. lower frequency for the x axis
% lowLim2  -double. lower frequency for the y axis
% q        -double. 1 for entropy & 2 for squared entropy
%% Outputs:
% Ent      -double. The calculated entropy
%-----------------------------------------------------------------------------------------------------------------
% Authors: Ioannis Ziogas & Charalampos Lamprou
% Copyright (C) 2022 Ioannis Ziogas and Charalampos Lamprou,SPBTU,ECE,AUTh
%-----------------------------------------------------------------------------------------------------------------

lowLim1 = ceil(lowLim1); lowLim2 = ceil(lowLim2); 
upLim1 = ceil(upLim1); upLim2 = ceil(upLim2);
%try
mat = abs(mat(lowLim1:upLim1,lowLim2:upLim2)).^q;
%catch
%    a = 1;
%end
p = mat./sum(sum(mat));
if isnan(p)
    p = zeros(size(p,1),size(p,2));   
end
p = p + 10e-10;
Ent = -sum(sum(p.*log(p)));

end