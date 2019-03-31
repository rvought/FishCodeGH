function [ output_args ] = all4thesize(idx, dat)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% size is idx(:,2) and stim is idx(:,1);
% PICK ALL STIM FOR A SIZE

availSizes = unique(idx(:,2));

fprintf('The available sizes are %i \n', availSizes);
sizeIwant = input('Pick a size (1-8): ');

stimIwant=idx(:,1);


% figure(1); clf;
%     ax(1) = subplot(212); plot(dat(sizeIwant).stim(stimIwant).ptim, dat(sizeIwant).stim(stimIwant).pos);

    ax(2) = subplot(211); hold on;
for pp = 1:length(dat(sizeIwant).stim(48).spiketimes) % Usually once, but sometimes more repetitions of same stimulus
    for k=1:length(dat(1).stim(48).spiketimes{1}); 
        plot(dat(1).stim(stimIwant).spiketimes{1}(k), pp-1, 'r*'); 
    end;
end;


end

