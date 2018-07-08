function out = cavefisheventdetector(data, in)
% Usage out = fisheventdtector(data, in)

pre = 200; post = 500;

numcolors = 50;
oclrs = lines(numcolors); % Colors for plotting up to numcolors different fishes

for j=1:length(data)
        
    z = zeros(1,length(data(j).tim));
    z(data(j).freq > (in(j).mmeanfreq + (2 * in(j).stdfreq))) = 1;
    
    out(j).eventidxs = find(diff(z) == 1);   
    
end

figure(27); clf; hold on;
for j = 1:length(data)
    
    if ~isempty(out(j).eventidxs)
    for k = 1:length(out(j).eventidxs)
        strt = max([out(j).eventidxs(k)-pre, 1]);
        ed = min([out(j).eventidxs(k)+post, length(data(j).tim)]);
       plot(data(j).tim(strt:ed), data(j).freq(strt:ed), 'Color', oclrs(k,:));
       
    end
    end
end



