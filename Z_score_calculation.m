% AptaZ algorithm part 1: Calculation of individual Z score
% Author: Daniel Wang
% Version: 1.0
% Updated: 2023-03-28

clc
clear all
close all

%% Parameter setting
psedocount = 5;

%% Make new folder for data storage
mkdir Z-results

%% Read and data preprocessing
% Read control sample
[file, path] = uigetfile('*.txt','Select the REFERENCE txt');
input_seq = readcell([path '\' file]);
[r c] = size(input_seq);
% Char conversion
for i = 1:length(input_seq)
    a = input_seq(i,c);
    c_seq(i,1) = string(a); %c_seq: detailed sequences in the control sample
end
% Num conversion
for i = 1:length(input_seq)
    a = regexp(input_seq(i,c-1),'\d*','match');
    a = a{1,1};
    b = a(1,length(a));
    c_count(i,1) = str2double(b); %c_counts: corresponding counts of individual sequences in the control sample
end

clear input_seq r c

% Read sorted samples
path = uigetdir(pwd,'Select the folder containing SORTED txt files');
file_list = dir(path);
j = length(file_list);
for k = 3:length(file_list)
    file = file_list(k).name;
    input_seq = readcell([path '\' file]);
    [r c] = size(input_seq);
    % Char conversion
    for i = 1:length(input_seq)
        a = input_seq(i,c);
        s_seq(i,1) = string(a); %s_seq: detailed sequences in the sorted sample
    end
    % Num conversion
    for i = 1:length(input_seq)
        a = regexp(input_seq(i,c-1),'\d*','match');
        a = a{1,1};
        b = a(1,length(a));
        s_count(i,1) = str2double(b); %s_counts: corresponding counts of individual sequences in the sorted sample
    end
    clear input_seq
    
    %% Calculate FC and Z score
    % Normalization to counts per million
    c_count_norm = c_count/sum(c_count)*1E6;
    s_count_norm = s_count/sum(s_count)*1E6;
    
    % Calculate Fold of change (FC)
    for i = 1:length(s_seq)
        [m,n] = find(c_seq == s_seq(i,1));
        if length(m) == 0
            norm_c(i) = psedocount; %if not match, assign the default psedocount
        else
            norm_c(i) = c_count_norm(m,1) + psedocount; %if matched, calculate the actual count
        end
        s_fc(i,1) = (s_count_norm(i) + psedocount)/norm_c(i); %calculate the fold change
    end
    % Calculate Z
    s_median = median(s_fc);
    s_fc_norm_log = log2(s_fc/s_median); %Z score
    
    % Rank the sequences based on Z score
    [s_fc_nls I] = sort(s_fc_norm_log,'descend');
    for i = 1:length(I)
        s_seq_s(i) = s_seq(I(i),1);
    end
    s_seq_s = s_seq_s';
    
    %% Save calculated matrix for Sum Z calculation
    i = length(file);
    disp(['completed ' num2str(k-2) ' file!']); 
    save([pwd '\Z-results\', file(1:i-3) 'mat'],'s_seq_s','s_fc_nls');
    
    %% Clean up
    clear s_seq norm_c s_fc
    clear s_median s_fc_norm_log
    clear s_fc_nls I
    clear s_seq_s
end
disp(['done']); 