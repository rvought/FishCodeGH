function [ strt, stp ] = brasilsamplelist(id, fish)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% CAVE
if id(1) == 1 % CAVE
   
    if id(2) == 3       
        strt = 0; stp = 1200;                
    end
    
    if id(2) == 4
          strt = 0; stp = 1200;        
        if sum(ismember(fish, 5)) > 0
          stp = 1000; 
        end
        if sum(ismember(fish, 6)) > 0
          strt = 400; 
        end
    end
    
    if id(2) == 5       
        strt = 0; stp = 1200;                
    end
    
    if id(2) == 6 
        strt = 0; stp = 1100;
        if sum(ismember(fish, 2)) > 0
          stp = 1000;
        end
        if sum(ismember(fish, 5)) > 0
          strt = 100;
        end
    end

    if id(2) == 7
        strt = 0; stp = 900;
        if sum(ismember(fish, 2)) > 0
          strt = 150; 
        end
        if sum(ismember(fish, 6)) > 0
          stp = 400;
        end
    end

    if id(2) == 8
        strt = 0; stp = 1250;
        if sum(ismember(fish, 2)) > 0
          strt = 1100; 
        end
    end    
    
    if id(2) == 10       
        strt = 0; stp = 900;                
    end
    
    if id(2) == 11       
        strt = 0; stp = 900;                
    end

    if id(2) == 12       
        strt = 0; stp = 1100;                
    end
    
    if id(2) == 13
        strt = 0; stp = 1200;
        if sum(ismember(fish, 8)) > 0
          strt = 300; 
        end
        if sum(ismember(fish, 9)) > 0
          stp = 300;
        end
        if sum(ismember(fish, 7)) > 0
          strt = 0; stp = 0;
        end
    end
    
    if id(2) == 14
        strt = 0; stp = 1000;
        if sum(ismember(fish, 7)) > 0
          strt = 300; 
        end
    end    
    

end
%% SURFACE
if id(1) == 2 % Surface

% Sample 1    
    if id(2) == 1
        strt = 0; stp = 600;
        
        if sum(ismember(fish, 3)) > 0
          stp = min([stp, 240]); 
        end
        
        if sum(ismember(fish, 4)) > 0
          strt = max([strt, 300]); 
          stp = min([stp, 500]); 
        end
        
        if sum(ismember(fish, 5)) > 0
          strt = max([strt, 550]); 
        end

        if sum(ismember(fish, 8)) > 0
          strt = max([strt, 100]); 
          stp = min([stp, 200]); 
        end
        
        if sum(ismember(fish, 10)) > 0
          strt = max([strt, 250]); 
          stp = min([stp, 550]); 
        end
        
    end
    
% Sample 2
    if id(2) == 2
        strt = 0; stp = 1050;
        
        if sum(ismember(fish, 1)) > 0
          strt = max([strt, 100]); 
          stp = min([stp, 500]); 
        end

        if sum(ismember(fish, 2)) > 0
          strt = max([strt, 700]); 
        end
    
        if sum(ismember(fish, 4)) > 0
          strt = max([strt, 350]); 
        end
        
        if sum(ismember(fish, 5)) > 0
          strt = max([strt, 350]); 
        end
        
        if sum(ismember(fish, 6)) > 0
          strt = max([strt, 300]); 
          stp = min([stp, 700]); 
        end
        
        if sum(ismember(fish, 7)) > 0
          strt = max([strt, 850]); 
        end
        
        if sum(ismember(fish, 8)) > 0
          strt = max([strt, 600]); 
          stp = min([stp, 600]); 
        end
        
        if sum(ismember(fish, 9)) > 0
          strt = max([strt, 600]); 
          stp = min([stp, 600]); 
        end
        
        if sum(ismember(fish, 10)) > 0
          strt = max([strt, 850]); 
        end
        
        if sum(ismember(fish, 12)) > 0
          strt = max([strt, 850]); 
        end
        
        if sum(ismember(fish, 14)) > 0
          stp = min([stp, 500]); 
        end

        if sum(ismember(fish, 15)) > 0
          strt = max([strt, 850]); 
        end
        
        if sum(ismember(fish, 16)) > 0
          strt = max([strt, 850]); 
        end
        
        if sum(ismember(fish, 17)) > 0
          strt = max([strt, 300]); 
          stp = min([stp, 700]); 
        end
        
        if sum(ismember(fish, 20)) > 0
          strt = max([strt, 600]); 
          stp = min([stp, 600]); 
        end
        if sum(ismember(fish, 21)) > 0
          strt = max([strt, 600]); 
          stp = min([stp, 600]); 
        end
    
    end        
        
% Sample 3
    if id(2) == 3
        strt = 0; stp = 600;

        if sum(ismember(fish, 2)) > 0
          stp = min([stp, 500]); 
        end
        
        if sum(ismember(fish, 4)) > 0
          strt = max([strt, 80]); 
          stp = min([stp, 180]); 
        end
        
        if sum(ismember(fish, 5)) > 0
          strt = max([strt, 200]); 
        end
        
        if sum(ismember(fish, 6)) > 0
          strt = max([strt, 600]); 
          stp = min([stp, 600]); 
        end
        
        if sum(ismember(fish, 7)) > 0
          stp = min([stp, 450]); 
        end
        
        if sum(ismember(fish, 8)) > 0
          strt = max([strt, 250]); 
        end

        if sum(ismember(fish, 9)) > 0
          strt = max([strt, 220]); 
        end
        
        if sum(ismember(fish, 10)) > 0
          strt = max([strt, 450]); 
        end
        
        if sum(ismember(fish, 11)) > 0
          strt = max([strt, 400]); 
        end
        
        if sum(ismember(fish, 12)) > 0
          stp = min([stp, 250]); 
        end
        
        if sum(ismember(fish, 13)) > 0
          strt = max([strt, 600]); 
          stp = min([stp, 600]); 
        end
        
        if sum(ismember(fish, 14)) > 0
          strt = max([strt, 600]); 
          stp = min([stp, 600]); 
        end
        
        if sum(ismember(fish, 15)) > 0
          strt = max([strt, 250]); 
          stp = min([stp, 320]); 
        end

        if sum(ismember(fish, 16)) > 0
          strt = max([strt, 50]); 
          stp = min([stp, 550]); 
        end
        
        if sum(ismember(fish, 17)) > 0
          strt = max([strt, 600]); 
          stp = min([stp, 600]); 
        end
        
        if sum(ismember(fish, 18)) > 0
          strt = max([strt, 600]); 
          stp = min([stp, 600]); 
        end
        
        if sum(ismember(fish, 19)) > 0
          strt = max([strt, 600]); 
          stp = min([stp, 600]); 
        end
        
        if sum(ismember(fish, 20)) > 0
          strt = max([strt, 440]); 
        end
        
        if sum(ismember(fish, 21)) > 0
          strt = max([strt, 540]); 
        end
        
        if sum(ismember(fish, 23)) > 0
          strt = max([strt, 600]); 
          stp = min([stp, 600]); 
        end
        
        if sum(ismember(fish, 24)) > 0
          strt = max([strt, 450]); 
        end
        
        if sum(ismember(fish, 25)) > 0
          strt = max([strt, 320]); 
          stp = min([stp, 550]); 
        end
        
        if sum(ismember(fish, 27)) > 0
          strt = max([strt, 250]); 
          stp = min([stp, 400]); 
        end
        
        if sum(ismember(fish, 29)) > 0
          stp = min([stp, 150]); 
        end
        
    end
    
% Sample 4
    if id(2) == 4
        strt = 0; stp = 620;

        if sum(ismember(fish, 1)) > 0
          strt = max([strt, 150]); 
          stp = min([stp, 250]); 
        end
        
        if sum(ismember(fish, 6)) > 0
          strt = max([strt, 470]); 
        end
        
        if sum(ismember(fish, 7)) > 0
          strt = max([strt, 380]); 
        end
       
        if sum(ismember(fish, 8)) > 0
          stp = min([stp, 160]); 
        end
        
        if sum(ismember(fish, 9)) > 0
          stp = min([stp, 80]); 
        end
        
        if sum(ismember(fish, 10)) > 0
          strt = max([strt, 220]); 
          stp = min([stp, 440]); 
        end
        
        if sum(ismember(fish, 13)) > 0
          strt = max([strt, 60]); 
          stp = min([stp, 220]); 
        end
        
        if sum(ismember(fish, 14)) > 0
          strt = max([strt, 440]); 
          stp = min([stp, 580]); 
        end
        
        if sum(ismember(fish, 15)) > 0
          strt = min([strt, 160]); 
        end
        
        if sum(ismember(fish, 17)) > 0
          strt = min([strt, 140]); 
        end

        if sum(ismember(fish, 18)) > 0
          strt = min([strt, 100]); 
        end
        
        if sum(ismember(fish, 19)) > 0
          strt = min([strt, 80]); 
        end
        
        if sum(ismember(fish, 21)) > 0
          strt = min([strt, 150]); 
        end 
        
    end
    
% Sample 5
    if id(2) == 5
        strt = 0; stp = 1100;
    
        if sum(ismember(fish, 1)) > 0
          strt = max([strt, 80]); 
          stp = min([stp, 230]); 
        end
        
        if sum(ismember(fish, 2)) > 0
          strt = max([strt, 150]); 
          stp = min([stp, 800]); 
        end
        
        if sum(ismember(fish, 3)) > 0
          strt = max([strt, 750]); 
        end        
        
        if sum(ismember(fish, 5)) > 0
          strt = max([strt, 80]); 
        end        
        
        if sum(ismember(fish, 7)) > 0
          stp = min([stp, 560]); 
        end
        
        if sum(ismember(fish, 9)) > 0
          strt = max([strt, 100]); 
        end        
        
        if sum(ismember(fish, 13)) > 0
          strt = max([strt, 250]); 
          stp = min([stp, 600]); 
        end

        if sum(ismember(fish, 16)) > 0
          strt = max([strt, 150]); 
          stp = min([stp, 800]); 
        end
        
        if sum(ismember(fish, 17)) > 0
          stp = min([stp, 820]); 
        end
        
        if sum(ismember(fish, 18)) > 0
          strt = max([strt, 150]); 
          stp = min([stp, 360]); 
        end

        % INDEXES 21 / 24 in orig are amazing less than 5Hz case
        
        
    end
    
end


end

