function out = movingconcat(numfiles)

    OrigVars = {};

for i = 1:numfiles

    clear V*;    
    uiopen('load');
    
    if i == 1;
        AllVars = who;
        for k = 1:length(AllVars);
            if AllVars{k}(1) == 'V';
                OrigVars{end+1} = AllVars{k};
            end;
        end;
        
        fprintf('Found %i variables \n', length(OrigVars));
        whos
       
        for j = 1:length(OrigVars)
 
        eval(OrigVars{j})
        NewVarName{j} = input('Give the new variable name. \n');
        str = ['out.' NewVarName{j} ' = ' OrigVars{j} ';'];
        eval(str);
        clear eval(OrigVars{j});
   
                
        end   
        out.I.values = out.I.values'; 
        out.T.values = out.T.values';
        out.P.values = out.P.values';
        out.ST.times = out.ST.times';
        out.K.times = out.K.times';
    end

    
    
    if i > 1;
        
    LastTime = length(out.T.values)*out.T.interval;
    
    for j = 1:length(OrigVars)
       
        str = [NewVarName{j} ' = ' OrigVars{j} ';'];
        eval(str);
        clear eval(OrigVars{j});
        who
        
    end;
    
    out.I.values = [out.I.values I.values']; 
    out.T.values = [out.T.values T.values'];
    out.P.values = [out.P.values P.values'];
    out.ST.times = [out.ST.times (ST.times' + LastTime + out.T.interval)];
    out.K.times = [out.K.times (K.times' + LastTime + out.T.interval)];
        
    end
    
end
    out.I.length = length(out.I.values);
        out.I.tim = out.I.interval:out.I.interval:out.I.length*out.I.interval;
    out.T.length = length(out.T.values);
        out.T.tim = out.T.interval:out.T.interval:out.T.length*out.T.interval;
    out.P.length = length(out.P.values);
        out.P.tim = out.P.interval:out.P.interval:out.P.length*out.P.interval;
    out.ST.length = length(out.ST.times);
    out.K.length = length(out.K.times);
    
