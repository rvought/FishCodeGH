% Plot spike triggered averages for error_pos, error_vel, error_acc, and error_jerk
% Load data first! Relies on iu_sta.m

% Load your data first (downsampled_data.mat)

% spks = spikes.times(spikes.times > 30 & spikes.times < 360);
% rspks = spikes_rand.times(spikes_rand.times > 30 & spikes_rand.times < 360);
% spks = spikes.times(spikes.times < 30);
% rspks = spikes_rand.times(spikes_rand.times < 30);
spks = spikes.times;
rspks = spikes_rand.times;

%% Calculate spike triggered averages
    fprintf('Calculating error_pos STA.\n');
    epos = iu_sta(spks, rspks, fish_pos, fs, 2);
    fprintf('Calculating error_vel STA.\n');
    evel = iu_sta(spks, rspks, fish_vel, fs, 2);
    fprintf('Calculating error_acc STA.\n');
    eacc = iu_sta(spks, rspks, fish_acc, fs, 2);
    fprintf('Calculating error_jerk STA.\n');
    ejerk = iu_sta(spks, rspks, fish_jerk, fs, 2);
    fprintf('And we are done!!!\n');

    %% Plot them all in one figure
    figure(1); clf; 

    subplot(2,2,1); title('Position'); hold on;
    plot([0, 0], [min(epos.MEAN), max(epos.MEAN)], 'k-', 'LineWidth',1);
    plot(epos.time, epos.MEAN, 'b-', 'LineWidth', 3);
    plot(epos.time, epos.randMEAN,'r-','LineWidth',3);

    subplot(222); title('Acceleration'); hold on;
    plot([0, 0], [min(eacc.MEAN), max(eacc.MEAN)], 'k-', 'LineWidth',1);
    plot(eacc.time, eacc.MEAN, 'b-', 'LineWidth', 3);
    plot(eacc.time, eacc.randMEAN,'r-','LineWidth',3);

    subplot(223); title('Velocity'); hold on;
    plot([0, 0], [min(evel.MEAN), max(evel.MEAN)], 'k-', 'LineWidth',1);
    plot(evel.time, evel.MEAN, 'b-', 'LineWidth', 3);
    plot(evel.time, evel.randMEAN,'r-','LineWidth',3);
    
    subplot(224); title('Jerk'); hold on;
    plot([0, 0], [min(ejerk.MEAN), max(ejerk.MEAN)], 'k-', 'LineWidth',1);
    plot(ejerk.time, ejerk.MEAN, 'b-', 'LineWidth', 3);
    plot(ejerk.time, ejerk.randMEAN,'r-','LineWidth',3);

%% Plot the error

figure(7); clf; 

    subplot(221); title('Position'); hold on;
    plot([0, 0], [min(epos.STD), max(epos.STD)], 'k-', 'LineWidth',1);
    plot(epos.time, epos.STD, 'b-', 'LineWidth', 3);
    plot(epos.time, epos.randSTD, 'r-', 'LineWidth', 3);

    subplot(222); title('Acceleration'); hold on;
    plot([0, 0], [min(eacc.STD), max(eacc.STD)], 'k-', 'LineWidth',1);
    plot(eacc.time, eacc.STD, 'b-', 'LineWidth', 3);
    plot(eacc.time, eacc.randSTD, 'r-', 'LineWidth', 3);

    subplot(223); title('Velocity'); hold on;
    plot([0, 0], [min(evel.STD), max(evel.STD)], 'k-', 'LineWidth',1);
    plot(evel.time, evel.STD, 'b-', 'LineWidth', 3);
    plot(evel.time, evel.randSTD, 'r-', 'LineWidth', 3);

    subplot(224); title('Jerk'); hold on;
    plot([0, 0], [min(ejerk.STD), max(ejerk.STD)], 'k-', 'LineWidth',1);
    plot(ejerk.time, ejerk.STD, 'b-', 'LineWidth', 3);
    plot(ejerk.time, ejerk.randSTD, 'r-', 'LineWidth', 3);

    