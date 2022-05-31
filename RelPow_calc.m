%% Thesis - Charalampos Lamprou 9114 & Ioannis Ziogas 9132 - AUTh ECE
% Thesis Project: Classification and Characterization of the Effect of Migraine 
% through Functional Connectivity Characteristics: Application to EEG 
% Recordings from a Multimorbid Clinical Sample

function V = RelPow_calc(cut_mat,total_mat,symmetric)

% This function calculates the relative power of a specific frequency block.
%% Inputs:
% cut_mat     -double matrix. Part of the whole matrix (total_mat)
% total_mat   -double matrix. The whole matrix
% symmetric   -double. 1 or 0. 1 if the total mat is symmetric with respect
%              to the x=y line and 0 in the opposite case.
%% Outputs:
% V           -double. The calculated relative power (as a percentage).

%-----------------------------------------------------------------------------------------------------------------
% Authors: Ioannis Ziogas & Charalampos Lamprou
% Copyright (C) 2022 Ioannis Ziogas and Charalampos Lamprou,SPBTU,ECE,AUTh
%-----------------------------------------------------------------------------------------------------------------

total_mat = abs(total_mat);
if symmetric
    total_mat = triu(total_mat);
    cut_mat = triu(abs(cut_mat));
end

Vtotal = trapz(trapz(total_mat.^2,2));
V = trapz(trapz(cut_mat.^2,2));
V = V/Vtotal*100;

end
