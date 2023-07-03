% AptaZ algorithm part 1: Calculation of individual Z score (fast approach)
% Author: Daniel Wang
% Version: 1.1
% Updated: 2023-07-03

clc
clear
close all

%% Parameter setting
psedocount = 5;

%% Make new folder for data storage
mkdir Z-results

%% Read and data preprocessing
% Read control sample
[file, path] = uigetfile('*.txt');
input_seq = readcell([path '\' file]);
[r, c] = size(input_seq);
% Char conversion
c_seq = string(zeros(length(input_seq),1));
for i = 1:length(input_seq)
    a = input_seq(i,c);
    c_seq(i,1) = string(a);
end
% Num conversion
c_count = zeros(length(input_seq),1);
for i = 1:length(input_seq)
    a = regexp(input_seq(i,c-1),'\d*','match');
    a = a{1,1};
    b = a(1,length(a));
    c_count(i,1) = str2double(b);
end

clear input_seq r c

% Read sorted samples
path = uigetdir;
file_list = dir(path);
j = length(file_list);
for k = 3:length(file_list)
    file = file_list(k).name;
    input_seq = readcell([path '\' file]);
    [r c] = size(input_seq);
    % Char conversion
    s_seq = string(zeros(length(input_seq),1));
    for i = 1:length(input_seq)
        a = input_seq(i,c);
        s_seq(i,1) = string(a);
    end
    % Num conversion
    s_count = zeros(length(input_seq),1);
    for i = 1:length(input_seq)
        a = regexp(input_seq(i,c-1),'\d*','match');
        a = a{1,1};
        b = a(1,length(a));
        s_count(i,1) = str2double(b);
    end
    clear input_seq
    
    %% Calculate FC and Z
    % Normalization to counts per million
    c_count_norm = c_count/sum(c_count)*1E6;
    s_count_norm = s_count/sum(s_count)*1E6;
    
    % Calculate Fold of change (FC)
    inter_seq = intersect(c_seq,s_seq);
    norm_c = zeros(length(s_seq),1);

    for j = 1:length(inter_seq)
        [m_c, ~] = find(c_seq == inter_seq(j,1),1);
        [m_s, ~] = find(s_seq == inter_seq(j,1),1);
        norm_c(m_s,1) = c_count_norm(m_c,1);
    end
    
    norm_c = norm_c + psedocount;
    s_fc = (s_count_norm + psedocount)./(norm_c);
    
    % Calculate Z
    s_median = median(s_fc);
    s_fc_norm_log = log2(s_fc/s_median);
    [s_fc_nls I] = sort(s_fc_norm_log,'descend');
    s_seq_s = string(zeros(length(I),1));
    for i = 1:length(I)
        s_seq_s(i) = s_seq(I(i),1);
    end
    
    %% Save for comprehensive analysis
    i = length(file);
    disp(['completed ' num2str(k-2) ' file!']); 
    save([pwd '\Z-results\', file(1:i-3) 'mat'],'s_seq_s','s_fc_nls');
    clear s_count c_count_norm s_count_norm
    clear inter_seq norm_c m_c m_s s_fc
    clear s_median s_fc_norm_log s_fc_nls I
    clear s_seq_s s_seq
end
disp(['done']); 