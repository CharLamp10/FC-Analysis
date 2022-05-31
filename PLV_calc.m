%% Thesis - Charalampos Lamprou 9114 & Ioannis Ziogas 9132 - AUTh ECE
% Thesis Project: Classification and Characterization of the Effect of Migraine 
% through Functional Connectivity Characteristics: Application to EEG 
% Recordings from a Multimorbid Clinical Sample

function plv_mean = PLV_calc(x, y, surr_struct)
% Compute the Phase Locking Value between two signals across trials, according to Lachaux, 
% Rodriguez, Martinerie, and Varela (1999). The PLV value ranges from 0, indicating random 
% phase differences, to 1 indicating a fixed phase difference. 
% phase_sig1 and phase_sig2 should be the phase values of the signals in radians, arranged as
% Samples x Trials. These can bed
% computed using the Wavelet or Hilbert transform, for example:
% phase_sig = angle(hilbert(BPS)); 
% Where BPS is the signal after band-pass filtering around the frequency range of interest. 
% 
%% Inputs:
% x             -double array. Signal
% y             -double array. Signal
% surr_struct   -struct. Contains info about the surrogate method. For more
%                info check ConnectivityAnalysis.m
%% Outputs:
% plv_mean      -double. The plv value between x and y.
% 
%-----------------------------------------------------------------------------------------------------------------
% Authors: Ioannis Ziogas & Charalampos Lamprou
% Copyright (C) 2022 Ioannis Ziogas and Charalampos Lamprou,SPBTU,ECE,AUTh
%-----------------------------------------------------------------------------------------------------------------

phase_x = angle(hilbert(x));
phase_y = angle(hilbert(y));
e = exp(1i*(phase_x - phase_y));
method = surr_struct.st_sig_method;
if method ~= "noSurrogate"
    alpha = surr_struct.alpha;
    N = surr_struct.N;
end
if method == "permutation"
    for i = 1:N
        ind = randperm(length(x));
        z = phase_y(ind);
        surr_e = exp(1i*(phase_x - z));
        plv_surr(i) = abs(mean(surr_e));
    end
    plv_thres = prctile(plv_surr,(1-alpha)*100);
    plv_mean = abs(mean(e));
    if plv_mean < plv_thres
        plv_mean = 0;
    end
elseif method == "block-resampling"
     blocksize = surr_struct.blocksize; % enter requred block size
     L_x = length(phase_x) - mod(length(phase_x),blocksize);  
     x_block = reshape(phase_x(1:L_x), blocksize, []); % divide the amp into blocks
     L_y = length(phase_y) - mod(length(phase_y),blocksize);  
     y_block = reshape(phase_y(1:L_y), blocksize, []); % divide the amp into blocks

     for i = 1:N   % calculated for N times
         random_xblock = randperm(size(x_block,1),1); %randomly select phase block number
         x_surr = x_block(random_xblock,:); %extract phase signal for random block
         random_yblock = randperm(size(y_block,1),1); %randomly select Amp block number
         y_ran = y_block(random_yblock,:); %extract Amp signal for random block
         y_surr = y_ran(randperm(length(y_ran)));
         e = exp(1i*(x_surr - y_surr));
         plv_surr(i) = e;
     end
    plv_thres = prctile(plv_surr,(1-alpha)*100);
    plv_mean = abs(mean(e));
    if plv_mean < plv_thres
        plv_mean = 0;
    end
elseif method == "noSurrogate"
    
    plv_mean = abs(mean(e));
end