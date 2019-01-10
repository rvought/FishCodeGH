%% DSI Calculations
function out = killerDSI(in)



%Direction
tail = 0;
head = 0;

for jj = 1:length(in.pos)
    if in.pos(jj) < -0.1
        tail = tail+1;
    end
    if in.pos(jj) > 0.1
        head = head + 1;
    end
end  

out.DSI = (head-tail)/(head+tail);
%out.DSI = (head-tail)/max([head tail]); %if we wanted to use the max
%instead of the sum
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

