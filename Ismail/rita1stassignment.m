%% 04_14_2019
load /Volumes/Cgate/Data2019/4_14_2019/downsampled_data.mat
j=1;
ismail(j).date=04142019;
ismail(j).filename='4_14_2019/downsampled_data.mat';

ismail(j).Fs = fs;
ismail(j).time = time;         

ismail(j).spikes = spikes;
ismail(j).spikes_rand = spikes_rand;

ismail(j).fish_pos = fish_pos;
ismail(j).fish_acc = fish_acc;
ismail(j).fish_vel = fish_vel;
ismail(j).fish_jerk = fish_jerk;

ismail(j).error_pos = error_pos;    
ismail(j).error_vel = error_vel;
ismail(j).error_acc = error_acc;
ismail(j).error_jerk = error_jerk;

ismail(j).shuttle_pos = shuttle_pos;
ismail(j).shuttle_vel = shuttle_vel;
ismail(j).shuttle_acc = shuttle_acc; 

%% 04_12_2019
j=2;
load /Volumes/Cgate/Data2019/4_14_2019/downsampled_data.mat

ismail(j).date=04122019;
ismail(j).filename='4_12_2019/brown2019_04_12_merged_wEric_ID123.mat';

ismail(j).Fs = fs;
ismail(j).time = time;         

ismail(j).spikes = spikes;
ismail(j).spikes_rand = spikes_rand;

ismail(j).fish_pos = fish_pos;
ismail(j).fish_acc = fish_acc;
ismail(j).fish_vel = fish_vel;
ismail(j).fish_jerk = fish_jerk;

ismail(j).error_pos = error_pos;    
ismail(j).error_vel = error_vel;
ismail(j).error_acc = error_acc;
ismail(j).error_jerk = error_jerk;

ismail(j).shuttle_pos = shuttle_pos;
ismail(j).shuttle_vel = shuttle_vel;
ismail(j).shuttle_acc = shuttle_acc; 
