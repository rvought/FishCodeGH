%% DSI Calculations
function out = objDSI(in)


%This codes only for seeing where the object spends the majority of its
%time

%Direction
lowrange = [];
highrange = [];

for jj = 1:length(in.objpos)
    if in.objpos(jj) < -0.1
        lowrange = [lowrange in.objpos(jj)];
    end
    if in.objpos(jj) > 0.1
        highrange = [highrange in.objpos(jj)];
    end
end  
dlow = sum(lowrange)/length(lowrange);
dhigh = sum(highrange)/length(highrange);
out.nDSI = (dhigh-abs(dlow))/(dhigh+abs(dlow));
fprintf('nDSI=')
fprintf(num2str(out.nDSI));
fprintf('\n')


%Velocity

lowvel = [];
highvel = [];

for kk = 1:length(in.objvel)
    if in.objvel(kk) < -0.002
        lowvel = [lowvel in.objvel(kk)];
    end
    if in.objvel(kk) > 0.002
        highvel = [highvel in.objvel(kk)];
    end
end  
vlow = sum(lowvel)/length(lowvel);
vhigh = sum(highvel)/length(highvel);
out.nVSI = (vhigh-abs(vlow))/(vhigh+abs(vlow));
fprintf('nVSI=')
fprintf(num2str(out.nVSI))
fprintf('\n')

%Acceleration
 
lowacc = [];
highacc = [];

for mm = 1:length(in.objacc)
    if in.objacc(mm) < -3E-6
        lowacc = [lowacc in.objacc(mm)];
    end
    if in.objacc(mm) > 3E-6
        highacc = [highacc in.objacc(mm)];
    end
end  
alow = sum(lowacc)/length(lowacc);
ahigh = sum(highacc)/length(highacc);
out.nASI = (ahigh-abs(alow))/(ahigh+abs(alow));
fprintf('nASI=')
fprintf(num2str(out.nASI))
fprintf('\n')
end

