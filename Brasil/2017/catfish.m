function out = catfish
% Use catfish to assemble grid data chunks

cntu = 0;
cntr = 0;

while cntu < 10

uiopen
uiopen

cntr = cntr+1;

if cntr == 1 % This is the first file - nothing needs to be done
    out = cleantime(tracks, particle);
    figure(1); clf; hold on; figure(2); clf; hold on;
    for j=1:length(out)
      figure(1); plot(out(j).tim, out(j).freq, '*');
      figure(2); plot(out(j).tx, out(j).ty, '*');
    end
end


if cntr > 1 % Concatenating a file
    
    % Get the last time point from the previously concatenated data
        for j=length(out):-1:1; maxtim(j) = out(j).tim(end); end;
        maxtim = max(maxtim);

    % Get the new data
        foo = cleantime(tracks, particle);
        
    % Now we have to align the data between recording frames
    
        if length(foo) == length(out) % Then we very likely have the same fish.

            figure(3); clf; hold on;
            for k = 1:length(out) % This is a pre-loop to compare fish frequencies
               prefreq = mean(out(k).freq); postfreq = mean(foo(k).freq); 
               prestd = std(out(k).freq); poststd = std(foo(k).freq);
               
               plot([1, 2], [prefreq, postfreq], '-*', 'LineWidth', 4);
               plot([1, 1], [prefreq-(5*prestd), prefreq+(5*prestd)], '-k', 'LineWidth', 2);
               plot([2, 2], [postfreq-(5*poststd), postfreq+(5*poststd)], '-k', 'LineWidth', 2);
               xlim([0.8 2.2]);   
            end

            % Do we have the same fish or not?
            
            % NEED TO FIX
            
            
        % Plot the concantenated data
        figure(1); clf; hold on; figure(2); clf; hold on;
            for j=1:length(out)
               out(j).freq = [out(j).freq foo(j).freq];
               out(j).tim = [out(j).tim foo(j).tim + maxtim];
               out(j).posidx = [out(j).posidx foo(j).posidx];
               out(j).tt = [out(j).tt foo(j).tt + maxtim];
               out(j).tx = [out(j).tx foo(j).tx];
               out(j).ty = [out(j).ty foo(j).ty];
                figure(1); plot(out(j).tim, out(j).freq, '*');
                figure(2); plot(out(j).tx, out(j).ty, '*');
            end
            
            % We will nevertheless check to match EOD frequencies
            % freq    posidx  tim     tt      tx      ty      x       y       

        end;
        
        if length(foo) ~= length(out) % We have a change in the number of fish. Joys.
            
            figure(1); clf; hold on; figure(2); clf; hold on;
            for j=1:length(out)
                figure(1); plot(out(j).tim, out(j).freq, '*');
                figure(2); plot(out(j).tx, out(j).ty, '*');
            end

            fixed = 0;
            
            figure(4); clf; 
            subplot(121); hold on;
            for j=1:length(out)
               plot(out(j).tim, out(j).freq, '*');
            end
            subplot(122); hold on;
            for j=1:length(foo)
               plot(foo(j).tim, foo(j).freq, '*');
            end            
            
            length(out)
            
            if length(foo) == length(out)-1 % We lost a fish.  This should be easier.
                
                lostfish = input('Which fish (starting from the top) did we lose? ');
                
                for j=1:lostfish-1
                out(j).freq = [out(j).freq foo(j).freq];
                out(j).tim = [out(j).tim foo(j).tim + maxtim];
                out(j).posidx = [out(j).posidx foo(j).posidx];
                out(j).tt = [out(j).tt foo(j).tt + maxtim];
                out(j).tx = [out(j).tx foo(j).tx];
                out(j).ty = [out(j).ty foo(j).ty];
                end;
                
                for j=lostfish+1:length(out)
                out(j).freq = [out(j).freq foo(j-1).freq];
                out(j).tim = [out(j).tim foo(j-1).tim + maxtim];
                out(j).posidx = [out(j).posidx foo(j-1).posidx];
                out(j).tt = [out(j).tt foo(j-1).tt + maxtim];
                out(j).tx = [out(j).tx foo(j-1).tx];
                out(j).ty = [out(j).ty foo(j-1).ty];
                end;
                
                fixed = 1 ;
            end
            
            if length(foo) == length(out)+1 % We gained a fish.  
                
                addfish = input('Which fish (starting from the top) did we add? ');
                length(out)
                out(end+1) = out(end);
                length(out)

                tmp = out;
                                
                if addfish < length(foo) % not the bottom fish
                    for j=addfish:length(foo)-1
                        out(j+1) = tmp(j);
                    end
                end
                
%                 if addfish == 1
%                     for j=2:length(out)
%                         tmp(j) = out(j);
%                     end
%                 end
                
                out(addfish).freq = [];
                out(addfish).tim = [];
                out(addfish).posidx = [];
                out(addfish).tt = [];
                out(addfish).tx = [];
                out(addfish).ty = [];
                
            for j=1:length(foo) % WILL THIS WORK???
               out(j).freq = [out(j).freq foo(j).freq];
               out(j).tim = [out(j).tim foo(j).tim + maxtim];
               out(j).posidx = [out(j).posidx foo(j).posidx];
               out(j).tt = [out(j).tt foo(j).tt + maxtim];
               out(j).tx = [out(j).tx foo(j).tx];
               out(j).ty = [out(j).ty foo(j).ty];
            end

                fixed = 1;
            end
            
            if fixed ~=1; fprintf('ARRRGHGHHHH!!!!! \n'); end;
            
            
        end;    


qwer = inputdlg('Done? If yes, 99 and then OK. If there are more files, just click OK:');

cntu = cntu + str2num(qwer{:});
if isempty(cntu) == 1; cntu = 0; end;

end


end




