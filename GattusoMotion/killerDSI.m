%% DSI Calculations
function out = killerDSI(in)



%Direction
lowrange = [];
highrange = [];

for jj = 1:length(in.pos)
    if in.pos(jj) < -0.1
        lowrange = [lowrange in.pos(jj)]; %should be 
    end
    if in.pos(jj) > 0.1
        highrange = [highrange in.pos(jj)];
    end
end  
dlow = sum(lowrange)/length(lowrange);
dhigh = sum(highrange)/length(highrange); %WHAT?
out.DSI = (dhigh-abs(dlow))/(dhigh+abs(dlow));
fprintf('DSI=')
fprintf(num2str(out.DSI));
fprintf('\n')


%Velocity

lowvel = [];
highvel = [];

for kk = 1:length(in.vel)
    if in.vel(kk) < -0.002
        lowvel = [lowvel in.vel(kk)];
    end
    if in.vel(kk) > 0.002
        highvel = [highvel in.vel(kk)];
    end
end  
vlow = sum(lowvel)/length(lowvel);
vhigh = sum(highvel)/length(highvel);
out.VSI = (vhigh-abs(vlow))/(vhigh+abs(vlow));
fprintf('VSI=')
fprintf(num2str(out.VSI))
fprintf('\n')

%Acceleration
 
lowacc = [];
highacc = [];

for mm = 1:length(in.acc)
    if in.acc(mm) < -3E-6
        lowacc = [lowacc in.acc(mm)];
    end
    if in.acc(mm) > 3E-6
        highacc = [highacc in.acc(mm)];
    end
end  
alow = sum(lowacc)/length(lowacc);
ahigh = sum(highacc)/length(highacc);
out.ASI = (ahigh-abs(alow))/(ahigh+abs(alow));
fprintf('ASI=')
fprintf(num2str(out.ASI))
fprintf('\n')
end

