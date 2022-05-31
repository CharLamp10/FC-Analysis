%% Thesis - Charalampos Lamprou 9114 & Ioannis Ziogas 9132 - AUTh ECE
% Thesis Project: Classification and Characterization of the Effect of Migraine 
% through Functional Connectivity Characteristics: Application to EEG 
% Recordings from a Multimorbid Clinical Sample

function [mean_feat, out_names] = pac_swarm_calc(param_struct)
% This routine is used for extraction of features from the PAC family.
% Various regional options,PAC measures,surrogate methods are included. 
% Two PAC calculation strategies are provided: 
%   A) Swarm Strategy: Swarm Decomposition derived Oscillations are treated
%   as separate signals, and PAC values are computed between low and high
%   frequency oscillations
%   B) Conventional Strategy: Swarm Decomposition Oscillations are summed
%   to form a signal that corresponds to a rhythm, thus resembling a filtered 
%   segmented signal. PAC comodulogram is computed between rhythm "filtered" signals
%   and original features for this comodulogram are calculated.
%   
%   AVAILABLE FEATURES FOR CONVENTIONAL STRATEGY: 
%           1) Relative Power -Calculate the PAC "power" percentage in each comodulogram
%           block, normalized by the whole PAC comodulogram "power"
%           2) Entropy1 - Calculate the Entropy value of each comodulogram block,
%           with reference the distribution of the specific block
%           3) TotalPAC - Calculate the Sum of all PAC values inside each
%           comodulogram block
%
%% Inputs
% param_struct   -contains all necessary initializations and parameters
%
%% Outputs
% mean_feat      - vector containing values of feature calculated
% out_names      - vector containing names of features calculated,
%                  corresponding 1-1 to mean_feat values

%-----------------------------------------------------------------------------------------------------------------
% Authors: Ioannis Ziogas & Charalampos Lamprou
% Copyright (C) 2022 Ioannis Ziogas and Charalampos Lamprou,SPBTU,ECE,AUTh
%-----------------------------------------------------------------------------------------------------------------


%% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
if param_struct.phAmpExtr == "tf"
    %Collect signals from input struct and sum all components of each electrode
    x = param_struct.a; %Region A electrodes
    y = param_struct.b; %Region B electrodes
    Fs = param_struct.Fs;
    %Low and High frequency intervals 
    param_struct.low = [param_struct.bandsLow(1),param_struct.bandsLow(end)];
    param_struct.high = [param_struct.bandsHigh(1),param_struct.bandsHigh(end)];

    if param_struct.regional == "intra-regional" || param_struct.regional == "anti-symmetric"
% Sum oscillations in each rhythm of each channel of each region x and y
% For each region, a cell with all oscillations of all electrodes summed is constructed
        x = gather_sigs_tf(x,[param_struct.bandsLow(1),param_struct.bandsLow(end)],Fs);
        y = gather_sigs_tf(y,[param_struct.bandsHigh(1),param_struct.bandsHigh(end)],Fs);
        param_struct.x_low_struct = x;
        param_struct.x_high_struct = y;
        [mean_feat,out_names] = calc_tf(param_struct);
    else
% Sum oscillations in each rhythm of each channel of each region x and y
% For each region, a cell with all oscillations of all electrodes summed is constructed
% Use each region once for low frequencies and once for high frequencies
        xl = gather_sigs_tf(x,[param_struct.bandsLow(1),param_struct.bandsLow(end)],Fs);
        yh = gather_sigs_tf(y,[param_struct.bandsHigh(1),param_struct.bandsHigh(end)],Fs);
        param_struct.x_low_struct = xl;
        param_struct.x_high_struct = yh;
        [mean_feat1,~] = calc_tf(param_struct);
        xh = gather_sigs_tf(x,[param_struct.bandsHigh(1),param_struct.bandsHigh(end)],Fs);
        yl = gather_sigs_tf(y,[param_struct.bandsLow(1),param_struct.bandsLow(end)],Fs);
        param_struct.x_low_struct = yl;
        param_struct.x_high_struct = xh;
        [mean_feat2,out_names] = calc_tf(param_struct);
        mean_feat = [mean_feat1,mean_feat2];
    end
else
    x = param_struct.a;
    y = param_struct.b;
    counter = 1;
    %Determine which rhythms will be used as "Low" Frequencies, by keeping
    %all other rhythms in "tempor" variable
    %Example: if wanted rhythms are 'delta','theta' -> param_struct.low = {'delta','theta'} ->
    % tempor = {'alpha','lo_beta','hi_beta','lo_gamma','hi_gamma'}
    % Then rhythms contained in tempor are deleted in gather_sigs_hil
    
    for t = 1:length(param_struct.rhythm_names)
        temporCounter = 0;
        for p = 1:length(param_struct.low)
            %Compare the cell of all rhythms, with the cell of desired
            %"Low" rhythms
            if ~strcmp(param_struct.rhythm_names{t},param_struct.low{p}) 
                temporCounter = temporCounter + 1;
            end
        end
        %If a rhythm does not match none of the "Low" rhythms, then add it
        %to tempor 
        if temporCounter == length(param_struct.low)
            tempor{counter} = param_struct.rhythm_names{t};
            counter = counter + 1;
        end
    end
    param_struct.Low = tempor;
    
    %"High" Rhythms
    counter = 1;
    for t = 1:length(param_struct.rhythm_names)
        temporCounter = 0;
        for p = 1:length(param_struct.high)
            %Compare the cell of all rhythms, with the cell of desired
            %"High" rhythms
            if ~strcmp(param_struct.rhythm_names{t},param_struct.high{p}) 
                temporCounter = temporCounter + 1;
            end
        end
        %If a rhythm does not match none of the "High" rhythms, then add it
        %to tempor (temporCounter w
        if temporCounter == length(param_struct.high)
            tempor{counter} = param_struct.rhythm_names{t};
            counter = counter + 1;
        end
    end
    param_struct.High = tempor;
    
    if param_struct.regional == "intra-regional" || param_struct.regional == "anti-symmetric"
% Classify oscillations in each rhythm of each channel of each region x and y
% For each region, a struct with all oscillations of all electrodes grouped is constructed
        xl = gather_sigs_hil(x,param_struct.Low,param_struct.rhythm_names);
        xh = gather_sigs_hil(y,param_struct.High,param_struct.rhythm_names);
        param_struct.xhigh = xh;
        param_struct.xlow = xl;
        [mean_feat,out_names] = calc_hil(param_struct);
    else
% Classify oscillations in each rhythm of each channel of each region x and y
% For each region, a struct with all oscillations of all electrodes grouped is constructed
% Use each region once for low frequencies and once for high frequencies
        xl = gather_sigs_hil(x,param_struct.Low,param_struct.rhythm_names);
        xh = gather_sigs_hil(y,param_struct.High,param_struct.rhythm_names);
        param_struct.xhigh = xh;
        param_struct.xlow = xl;
        [mean_feat1,~] = calc_hil(param_struct);
        
        xl = gather_sigs_hil(y,param_struct.Low,param_struct.rhythm_names);
        xh = gather_sigs_hil(x,param_struct.High,param_struct.rhythm_names);
        param_struct.xhigh = xh;
        param_struct.xlow = xl;
        [mean_feat2,out_names] = calc_hil(param_struct);
        
        mean_feat = [mean_feat1, mean_feat2];
    end
end
    

function [reg_table] = gather_sigs_tf(x,bandsLims,Fs)
% This function receives a region containing some electrodes as an input.
% Each electrode has some rhythms in form of oscillations. Oscillations that 
% belong in the specified band limits (bandsLims) are kept. Then these oscillations 
% are summed into one signal. This procedure is repeated channel-wise. If a channel
% has no oscillations that belong in the specified frequency interval, then
% it is filled with NaNs.
%
% INPUTS: x - struct, contains channels. Channels (structs) contain cells of oscillations
%         bandsLims - 1x2 array, limits of frequency interval [low up]
%         Fs - sampling rate
% OUTPUTS: reg_table - cell, contains summed oscillations for each rhythm
%

%% Written by: Ioannis Ziogas && Charalampos Lamprou, December 2021
    chans = fieldnames(x);
    count = 1;   
    %For each channel:
    for k = 1:length(chans)  
        %Oscillations belonging to the same rhythm are separated - unwrapped
        sig = unwrap_SwDs(x.(chans{k})); 
        sig1 = {};
        for s = 1:length(sig)
            flag = belongs(sig{s},bandsLims,Fs); %Find if signal belongs into specified bands limits
            if flag
                sig1{s} = sig{s}; 
            end
        end
%        sig1 = cell2table(sig1);
%         idx = all(cellfun(@isempty,sig{:,:}),2);
%         sig(idx,:)=[];
%        sig1 = table2array(sig1);
        if ~isempty(sig1)
            %If a channel has oscillations in this frequency interval
            sig1 = horzcat(sig1{:}); %Cell2array and remove empty cells
            [lx,hx] = size(sig1);
            if hx > lx
                sig1 = transpose(sig1);
                warning("Rhythms should be in shape [samples,realizations]")
            end                              
            reg_table(:,count) = sum(sig1,2);  
%             reg_chan_names(count) = chans{k};    
            count = count + 1;                       
        else
%             If a channel has 0 oscillations in this frequency interval
            temp = 0;
            for n = 1:length(sig)
                if ~isempty(sig{n})
                    temp = temp + sig{n};
                end
            end
            %Fill reg_table with NaNs
            reg_table(:,count) = NaN(length(temp),1);
            count = count + 1;
        end
    end
    %Then transform reg_table from array to cell
    %Final cell array contains summed oscillations for each rhythm
    if isstruct(reg_table)
        reg_table = struct2cell(reg_table);
        reg_table = unwrap_SwDs(reg_table); %if there are two SwDs in a band then separate them
    elseif ismatrix(reg_table)
        [l,h] = size(reg_table);
        if l < h
            reg_table = transpose(reg_table);
        end
        for jj= 1:length(reg_table(1,:))
            new_x{jj,1} = reg_table(:,jj);
        end
        reg_table = new_x;
    end
    
function flag = belongs(signal,limits,Fs)
% This function receives as inputs a signal and a bandwidth, and determines
% if the spectral content of the signal belongs in this bandwidth. 
% INPUTS: signal - numeric vector
%         limits - frequency limits of bandwidth, 1x2 array, [low high]
%         Fs     - sampling rate
%
% OUTPUTS: flag - binary flag, 1 if belongs, 0 otherwise

%% Written by: Ioannis Ziogas && Charalampos Lamprou, December 2021

        if ~isempty(signal)
            L = length(signal);
            nfft = 2^nextpow2(L); %Find optimal nfft 
            step = Fs/nfft; %Frequency domain step
            X = (abs(fft(signal,nfft))/L).^2; % FFT Power Spectrum
            X = X(1:end/2); % Discard symmetric side
            X_norm = X./sum(X); %Normalize the spectrum
            % Find spectral content percentage
            per = sum(X_norm(floor(limits(1)/step):floor(limits(2)/step)));
            if per > 0.8 
% 80 % of the signal's spectral content must be contained in the specified bandwidth
                flag = 1;
            else
                flag = 0;
            end  
        else
            flag = 0;
        end
end %End of belongs function
end %End of gather_sigs_tf function

function [mean_feat,out_names] = calc_tf(param_struct)
% This function receives as inputs two regions in the form of cells,
% containing oscillations of all region electrodes, grouped by rhythms.
% Then for each rhythm pair, the PAC comodulogram is computed, segmented into
% blocks indicated by specified bands. Features are calculated for each of
% these blocks. Three different PAC measures are available: MVL, MI, GLM
% Various surrogate data methods to ensure PAC values statistical
% significance are included.
%  
%INPUTS: param_struct - contains all necessary parameters and initializations
%
% OUTPUTS: mean_feat - vector containing all features calculated for all
%                   bands combinations
%          out_names - vector containing names of all features calculated.
%          1-1 correspondance to mean_feat vector

%% Written by: Ioannis Ziogas && Charalampos Lamprou, December 2021 


%% Initializations
    symmetric = 0; %This argument specifies if features like Relative Power 
    %are calculated over the whole comodulogram plane or only for half the
    %comodulogram plane
    type = param_struct.regional;
    x_low_cell = param_struct.x_low_struct; %Cell containing low rhythms of Region A 
    x_high_cell = param_struct.x_high_struct; %Cell containing low rhythms of Region B
    bandsLow = param_struct.bandsLow;
    bandsHigh = param_struct.bandsHigh;
    band_method = param_struct.band_method; %Conventional or non-conventional
    if band_method == "conventional"
        conv_bandsLow = param_struct.conv_bandsLow;
        conv_bandsHigh = param_struct.conv_bandsHigh;
    end
    features = param_struct.features; % Comodulogram Features to be calculated
    %Iterate through low frequency and high frequency signals
    for i1 = 1:length(x_low_cell)
        for j1 = 1:length(x_high_cell)
            if type == "single" || type == "intra-regional" || type == "anti-symmetric"
                ifcond = (j1 ~= i1);
            elseif type == "inter-regional"  || type == "left-right"
                ifcond = 1;
            end
            if ifcond && ~isempty(x_low_cell{i1}) && ~isempty(x_high_cell{j1}) && ~(all(isnan(x_low_cell{i1})) || all(isnan(x_high_cell{j1})))% j > i to discard duplicates
            % j1~=i1: x_low_cell and x_high_cell will be identical in these cases, so we
            % don't want to calculate pac between the same signal
                param_struct.x_low = x_low_cell{i1};
                param_struct.x_high = x_high_cell{j1};
                [pac_measure,~] = tf_PAC_surr(param_struct); %Returns comodulogram and max PAC value
                %Segment the comodulogram plane into frequency blocks
                for b1 = 1:2:length(bandsLow) - 1
                    %Find frequencies that correspond to f2 = [bands(i),bands(i+1)]
                    indLow = (b1+1)/2;
                    if band_method == "conventional"
                        indnLow = conv_bandsLow(indLow);
                    else
                        indnLow = string(indLow);
                    end
                    lowLim1 = (bandsLow(b1) - bandsLow(1))/param_struct.step_low + 1;
                    upLim1 = (bandsLow(b1+1) - bandsLow(1))/param_struct.step_low + 1;
                    pacTemp1 = zeros(size(pac_measure,1),size(pac_measure,2));
                    pacTemp1(:,lowLim1:upLim1) = pac_measure(:,lowLim1:upLim1);
                    %figure('Visible','off')
                    %contour(waxis,waxis,abs(bspecTemp1),4)

                    for b2 = 1:2:length(bandsHigh) - 1
                        %Find frequencies that correspond to f1 = [bands(j),bands(j+1)]
                        indHigh = (b2+1)/2;
                        if band_method == "conventional"
                            indnHigh = conv_bandsHigh(indHigh);
                        else
                            indnHigh = string(indHigh);
                        end
                        lowLim2 = (bandsHigh(b2) - bandsHigh(1))/param_struct.step_high + 1;
                        upLim2 = (bandsHigh(b2+1) - bandsHigh(1))/param_struct.step_high + 1;
                        pacTemp2 = zeros(size(pac_measure,1),size(pac_measure,2));
                        %pacTemp2 now is the region that corresponds to bandX -bandY
                        %coupling
                        pacTemp2(lowLim2:upLim2,:) = pacTemp1(lowLim2:upLim2,:);
                        %figure('Visible','off')
                        %contour(waxis,waxis,abs(bspecTemp2),4)

                        if sum(contains(features,"RelPow")) > 0
                    %If it is the first loop iteration, create RelPow variable, and create new struct field 
                            if ~exist('RelPow','var')
                                RelPow.(join(["RelPow_band",indnLow,"_",indnHigh],'')) = RelPow_calc(pacTemp2,pac_measure,symmetric);% RelPow_calc(pacTemp2,tf_MVL_all,symmetric);
                            elseif exist('RelPow','var') && isfield(RelPow,join(["RelPow_band",indnLow,"_",indnHigh],''))
                    %If it is not the first loop iteration, and current rhythm combination already exists 
                    %as a field, append RelPow value to already existing RelPow struct field
                                RelPow.(join(["RelPow_band",indnLow,"_",indnHigh],'')) = RelPow_calc(pacTemp2,pac_measure,symmetric);% RelPow_calc(pacTemp2,tf_MVL_all,symmetric);
                            else
                    %If it is not the first loop iteration, and current rhythm combination does not exist 
                    %as a field, create a new struct field in already existing RelPow struct 
                                RelPow.(join(["RelPow_band",indnLow,"_",indnHigh],'')) = RelPow_calc(pacTemp2,pac_measure,symmetric);% RelPow_calc(pacTemp2,tf_MVL_all,symmetric);
                            end
                        else 
                            %warning("RelativePower feature was not calculated. You should specify a feature name that contains RelPow if you want to calculate it")
                        end

                        if sum(contains(features,"Ent1")) > 0
                            %If it is the first loop iteration, create Ent1 variable, and create new struct field
                            if ~exist('Ent1','var')                                    
                                Ent1.(join(["Ent1_band",indnLow,"_",indnHigh],'')) = bspecEntropy(pacTemp2,upLim2,upLim1,lowLim2,lowLim1,1);                                   
                            elseif exist('Ent1','var') && isfield(Ent1,join(["Ent1_band",indnLow,"_",indnHigh],''))                      
                                %If it is not the first loop iteration, and current rhythm combination already exists 
                                %as a field, append Ent1 value to already existing Ent1 struct field    
                                Ent1.(join(["Ent1_band",indnLow,"_",indnHigh],''))(end+1) = bspecEntropy(pacTemp2,upLim2,upLim1,lowLim2,lowLim1,1);                            
                            else
                                %If it is not the first loop iteration, and current rhythm combination does not exist 
                                %as a field, create a new struct field in already existing Ent1 struct 
                                Ent1.(join(["Ent1_band",indnLow,"_",indnHigh],'')) = bspecEntropy(pacTemp2,upLim2,upLim1,lowLim2,lowLim1,1);                                    
                            end
                        else 
                            %warning("Entropy1 feature was not calculated. You should specify a feature name that contains Ent1 if you want to calculate it")
                        end

                        if sum(contains(features,"TotalPAC")) > 0
                            %If it is the first loop iteration, create TotalPAC variable, and create new struct field
                            if ~exist('TotalPAC','var')      
                                TotalPAC.(join(["TotalPAC_band",indnLow,"_",indnHigh],'')) = sum(sum(pacTemp2));                                   
                            elseif exist('TotalPAC','var') && isfield(TotalPAC,join(["TotalPAC_band",indnLow,"_",indnHigh],''))                      
                                %If it is not the first loop iteration, and current rhythm combination already exists 
                                %as a field, append TotalPAC value to already existing TotalPAC struct field    
                                TotalPAC.(join(["TotalPAC_band",indnLow,"_",indnHigh],''))(end+1) = sum(sum(pacTemp2));                                
                            else
                                %If it is not the first loop iteration, and current rhythm combination does not exist 
                                %as a field, create a new struct field in already existing TotalPAC struct
                                TotalPAC.(join(["TotalPAC_band",indnLow,"_",indnHigh],'')) = sum(sum(pacTemp2));                                     
                            end
                        else 
                            warning("TotalPAC feature was not calculated. You should specify a feature name that contains TotalPAC if you want to calculate it")
                        end
                    end
                end
            end
        end
    end
    out_names = [];
    %Initialize out_names vector
    if exist('RelPow','var') && ~isempty(RelPow)
        %If Relative Power was computed 
        names = fieldnames(RelPow);
        out_names = [out_names;names];  
        for i = 1:length(names)
            name = names{i};
            if sum(isnan(mean(RelPow.(name)))) >= 1
                error("RelPow contains NaN")
            else
                RelPow.(name) = mean(RelPow.(name));
            end
            mean_RelPow(i) = RelPow.(name);
        end
    elseif ~exist('RelPow','var')
        mean_RelPow = [];
    end

    if exist('Ent1','var') && ~isempty(Ent1)
        %If Entropy1 was computed 
        names = fieldnames(Ent1);
        out_names = [out_names;names];  
        for i = 1:length(names)
            name = names{i};
            if sum(isnan(mean(Ent1.(name)))) >= 1
                error("Ent1 contains NaN")
            else
                Ent1.(name) = mean(Ent1.(name));
            end
            mean_Ent1(i) = Ent1.(name);
        end
    elseif ~exist('Ent1','var')
        mean_Ent1 = [];
    end
    
    if exist('TotalPAC','var') && ~isempty(TotalPAC)
        %If Entropy1 was computed 
        names = fieldnames(TotalPAC);
        out_names = [out_names;names];  
        for i = 1:length(names)
            name = names{i};
            if sum(isnan(mean(TotalPAC.(name)))) >= 1
                error("TotalPAC contains NaN")
            else
                TotalPAC.(name) = mean(TotalPAC.(name));
            end
            mean_TotalPAC(i) = TotalPAC.(name);
        end
    elseif ~exist('TotalPAC','var')
        mean_TotalPAC = [];
    end
    
    if exist('mean_RelPow','var') || exist('mean_Ent1','var') || exist('mean_TotalPAC','var') 
        mean_feat = [mean_RelPow, mean_Ent1, mean_TotalPAC];
    end
    
    if isempty(mean_feat)
        mean_feat = NaN(1,param_struct.bands_length*length(param_struct.features));
    end
    
end %End of calc_tf function


%================================================================================


function [x] = gather_sigs_hil(x,HighLow, rhythms_names)
% This function receives a region struct x containing all regional electrodes.
% Each electrode is a cell containing oscillations grouped by rhythms.
% This function classifies oscillations inside each electrode into rhythm named fields.
% 
% INPUTS: x - struct, contains cell vectors - electrodes, that contain
%                       oscillations
%         HighLow - rhythm names that are considered for this region. Can
%                    be low frequency rhythms or high frequency rhythms (for PAC calculation)
%         rhythms_names - cell vector containing names of all possible
%                           rhythms
%
% OUTPUTS: x - struct, contains structs - electrodes, that contain
%               classified oscillations

%% Written by: Ioannis Ziogas && Charalampos Lamprou, December 2021

    chans = fieldnames(x);       
    for k = 1:length(chans)  
       if chans{k} ~= "NoChannels"
           x.(chans{k}) = cell2struct(x.(chans{k}), rhythms_names);
       end
    end
    for k = 1:length(chans) 
       if chans{k} ~= "NoChannels"
           x.(chans{k}) = rmfield(x.(chans{k}),HighLow);
       end
    end
     
end


function [mean_PAC,names] = calc_hil(par_struct)
%% Uses the hil_PAC_surr.m routine to calculate PAC value. In this approach, 
% Swarm Decomposition acquired oscillations are used and paired to obtain
% PAC values. Two regions are passed in the form of 'xhigh' and 'xlow'
% variables. The first region (xhigh) contains only high frequency oscillations
% and the second region (xlow) contains only low frequency oscillations,
% grouped in rhythms based on their spectral content.
% This function goes on and calculates PAC between low and high frequency
% oscillations, storing values to a 'mPAC' struct, in the form of vectors, one for each
% low-high oscillation pair - e.g. mPAC.delta_lo_gamma 
% The average PAC value of this pair is returned, and is boosted according to 
% a Boost Criterion. Names of all PAC pairs are also returned.
% Note: Regional Cases are handled differently
%
% INPUTS: par_struct - contains:
%                      xhigh & xlow: Region A and Region B electrodes, 
%                      pacThresh: significance Threshold (here set to 0)
%                      regional: Indicates regional case
%                      high & low: High and Low frequency rhythm names        
%                      bands_length: Number of rhythm combinations
%
% OUTPUTS: mean_PAC - 1xN vector: mean PAC value of each
%                   rhythm combination (of the N total combinations)
%                   for the A-B region pair
%          names - Nx1 cell vector: names of rhythm combinations (N total combs)
%% Written by: Charalampos Lamprou && Ioannis Ziogas , December 2021


    xhigh = par_struct.xhigh; %Fetch regional electrodes from which "High" rhythms will be used
    xlow = par_struct.xlow; %Fetch regional electrodes from which "Low" rhythms will be used
    electrodesHigh = fieldnames(xhigh); %Fetch "High" electrode names
    electrodesLow = fieldnames(xlow); %Fetch "Low" electrode names
    pacThresh = par_struct.pacThresh; %A PAC significance threshold - 
    %Here set to 0, adding importance to PAC values "passing" the surrogate test 

    elec_case = par_struct.regional; %Fetch regional case - e.g. inter-regional,intra-regional,...

%The "NoChannels" case occurs only when a region has 0 active electrodes
%If this case does not hold, then the procedure is executed: 
    if ~strcmp(electrodesHigh{1,1}, "NoChannels") && ~strcmp(electrodesLow{1,1}, "NoChannels")
        for eh = 1:length(electrodesHigh) %For each electrode 'el':
            %Take "High" rhythms of electrode 'el'
            xxhigh = xhigh.(electrodesHigh{eh}); %xx-electrode: struct, each field is a rhythm
            for el = 1:length(electrodesLow) %For each electrode of the electrodes that contain "Low" rhythms:
                if (elec_case == "intra-regional" || elec_case == "anti-symmetric" ...
                        || elec_case == "single")
                    ifcond = (el ~= eh);% el ~= eh: Because the "intra-regional" and "anti-symmetric" cases
%involve a single region, they have the same electrodes, and we don't want to measure PAC
%inside an electrode (between "Low" and "High" rhythms of the same electrode
                elseif elec_case == "inter-regional" || elec_case == "left-right"
                    ifcond = 1;
                end
                
                if ifcond
                    xxlow = xlow.(electrodesLow{el}); %xx-electrode: struct, each field is a rhythm
                    %Then iterate through chosen "High" and "Low" rhythms
                    for rh = 1:length(par_struct.high)  
                        for rl = 1:length(par_struct.low)
                            %xxx-rhythm: vector - can contain one oscillation,
                            %           or matrix - many oscillations 
                            xxxhigh = xxhigh.(par_struct.high{rh});                                
                            xxxlow = xxlow.(par_struct.low{rl});
                            %If the "High" and "Low" rhythms are not empty
                            if ~isempty(xxxhigh) && ~isempty(xxxlow)
                                %Count how many oscillations are present
                                lenh = size(xxxhigh,2); 
                                lenl = size(xxxlow,2);
                                %For each combination of rhythms, compute PAC
                                for lh = 1:lenh 
                                    for ll = 1:lenl
                                        if ~exist('mPAC','var')
                        %If it is the first loop iteration, create mPAC variable, and create new struct field 
                                            par_struct.x_low = xxxlow(:,ll); par_struct.x_high = xxxhigh(:,lh);
                                            [mPAC.(join(["mPAC",par_struct.low{rl},'_',par_struct.high{rh}],'')),~] ...    
                                                = hil_PAC_surr(par_struct); % calculate PAC                                   
                                        elseif exist('mPAC','var') && isfield(mPAC,join(["mPAC",par_struct.low{rl},'_',par_struct.high{rh}],''))
                        %If it is not the first loop iteration, and current rhythm combination already exists 
                        %as a field, append PAC value to already existing mPAC struct field                   
                                            par_struct.x_low = xxxlow(:,ll); par_struct.x_high = xxxhigh(:,lh);
                                            [mPAC.(join(["mPAC",par_struct.low{rl},'_',par_struct.high{rh}],''))(end+1), ~] ...
                                                = hil_PAC_surr(par_struct); % calculate PAC
                                        else  
                        %If it is not the first loop iteration, and current rhythm combination does not exist 
                        %as a field, create a new struct field in already existing mPAC struct 
                                            par_struct.x_low = xxxlow(:,ll); par_struct.x_high = xxxhigh(:,lh);
                                            [mPAC.(join(["mPAC",par_struct.low{rl},'_',par_struct.high{rh}],'')), ~] ...    
                                                = hil_PAC_surr(par_struct); % calculate PAC
                                        end  
                                    end
                                end
                            else
                        %If the "High" and/or "Low" rhythms are empty, then follow the previous
                        %procedure (described in the 382-399 'if' block) by filling struct fields with NaN's
                                if ~exist('mPAC','var')
                                    mPAC.(join(["mPAC",par_struct.low{rl},'_',par_struct.high{rh}],'')) ...    
                                        = NaN; % calculate PAC                                   
                                elseif exist('mPAC','var') && isfield(mPAC,join(["mPAC",par_struct.low{rl},'_',par_struct.high{rh}],''))
                                    mPAC.(join(["mPAC",par_struct.low{rl},'_',par_struct.high{rh}],''))(end+1) ...
                                        = NaN; % calculate PAC
                                else  
                                    mPAC.(join(["mPAC",par_struct.low{rl},'_',par_struct.high{rh}],'')) ...    
                                        = NaN; % calculate PAC
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    mean_PAC = NaN(1,par_struct.bands_length); %Initialize output of function
    if exist('mPAC','var')
%If mPAC variable exists (If PAC extraction procedure was executed successfully) 
        names = fieldnames(mPAC); %Fetch rhythm combinations calculated
        %For each combination:
        for i = 1:length(names)
            name = names{i};
%Check if all values in the rhythm combination vector are NaNs
            if all(isnan(mPAC.(name))) == 1
%If true, give this combination a mean PAC value of 0 indicating no coupling 
%between this particular pair of "Low" and "High" rhythms. 
%This case can only happen when either rhythm of the pair is empty
                mPAC.(name) = 0;
            else
%If both rhythms had at least one oscillation, find PAC values that pass the set significance threshold
                inds = find(mPAC.(name) > pacThresh); 
%% Boost Criterion
%Apply an Additive Boosting Factor to the mean PAC value of a combination
%Combinations that have A) many significant and B) strong couplings, are favored by this boost
%Combinations that have many significant couplings are favored by the first factor of the Boost
%Combinations that have strong couplings are favored by the second factor of the Boost
                pacBoost = (length(inds) / length(mPAC.(name))) *  mean(mPAC.(name)(inds),'omitnan');             
                if isempty(inds)
                %If no values pass the surrogate test, set mean PAC to 0
                    mPAC.(name) = 0;
                else
                %Apply the Boost to the mean PAC value of the
                %combination
                    mPAC.(name) = mean(mPAC.(name)(inds),'omitnan') + pacBoost;
                end
            end
            mean_PAC(i) = mPAC.(name);
        end
    elseif ~exist('mPAC','var') 
        %If no PAC was calculated, return an empty name vector and NaNs
        names = [];
        mean_PAC = NaN(1,par_struct.bands_length);
    end 
end %End of calc_hil function

end % End of pac_swarm_calc function
