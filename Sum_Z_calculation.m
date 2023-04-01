% AptaZ algorithm part 2: Calculation of Sum Z score
% Author: Daniel Wang
% Version: 1.0
% Updated: 2023-03-28

clc
clear all
close all

%% Make lists of parameters
flow_rate = zeros(1,1);
concent = zeros(1,1);
zone = zeros(1,1);
z_score = zeros(1,1);
z_seq = string(zeros(0,0));

%% Read samples
path = uigetdir(pwd,'Select the folder containing calculated Z scores');
file_list = dir(path);

%% Make new folder for data storage
mkdir Sum-Z-results

for k = 3:length(file_list)
    file = file_list(k).name;
    %% Read parameters and calculate the weights
    weight_temp = regexp(file,'\d*','match');
    flow_rate = str2double(weight_temp{1,1}); %weight by flow rate
    concent = str2double(weight_temp{1,2}); %weight by target concentration
    zone = str2double(weight_temp{1,3}); %weight by zone index
    weight = (flow_rate/16)*(9/(3^(zone-1)))*(6-log10(concent)); %overall weight
    
    %% Read seq and z
    load([path '\' file]);
    % Format check and correction
    if size(s_fc_nls,1) < size(s_fc_nls,2)
        s_fc_nls = s_fc_nls';
    end
    if size(s_seq_s,1) < size(s_seq_s,2)
        s_seq_s = s_seq_s';
    end
    
    %% Calculate Sum Z score
    [exist_boolean, exist_index] = ismember(s_seq_s,z_seq);
    j = length(z_seq);
    for i = 1:length(s_fc_nls)
        if exist_boolean(i) == 0 %if the sequence is new, add it to the pool of calculation
            j = j + 1;
            z_seq(j,1) = s_seq_s(i,1); %z_seq: all sequences presented in any of the sorted samples
            z_score(j,1) = weight*s_fc_nls(i,1); %z_score: corresponding Sum Z score for individual sequences
        else %if the sequence is not new, add the individual Z score to the existing Sum Z score
            z_score(exist_index(i),1) = z_score(exist_index(i),1) + weight*s_fc_nls(i,1);
        end
    end
    clear s_seq_s
    clear s_fc_nls
    disp(['completed ' num2str(k-2) ' file!']); 
end

%% Rank sequence based on Sum Z score
[z_score_rank, z_score_index] = sort(z_score,'descend');
z_seq_rank = string(zeros(0,0));
for i = 1:length(z_score_index)
    z_seq_rank(i,1) = z_seq(z_score_index(i));
end
figure
plot(z_score_rank,'LineWidth',10);
xlabel('Index of sequence');
ylabel('Z score');

writematrix(z_seq_rank,[pwd '\Sum-Z-results\' 'z_seq_rank.csv']);
writematrix(z_score_rank,[pwd '\Sum-Z-results\' 'z_score_rank.csv']);

disp('done');