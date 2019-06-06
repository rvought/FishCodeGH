function [ strt, stp ] = brasilsamplelist(id, fish)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


strt = []; stp = [];

if id(1) == 1 % CAVE
   
    if id(2) == 3       
        strt = 0; stp = 1200;                
    end
    
    if id(2) == 4
        if sum(ismember(fish, 6)) > 0
          strt = 400; stp = 1200;
        else
          strt = 0; stp = 1200;
        end
    end
    
    if id(2) == 5       
        strt = 0; stp = 1200;                
    end
    
    if id(2) == 6
        if sum(ismember(fish, 2)) > 0
          strt = 0; stp = 1000;
        elseif sum(ismember(fish, 5)) > 0
          strt = 100; stp = 1100;
        end
    end

    if id(2) == 7
        if sum(ismember(fish, 2)) > 0
          strt = 150; stp = 900;
        elseif sum(ismember(fish, 6)) > 0
          strt = 0; stp = 400;
        else
          strt = 0; stp = 900;
        end
    end

    if id(2) == 8
        if sum(ismember(fish, 2)) > 0
          strt = 1100; stp = 1250;
        else
          strt = 0; stp = 1250;
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
        if sum(ismember(fish, 8)) > 0
          strt = 300; stp = 1200;
        elseif sum(ismember(fish, 9)) > 0
          strt = 0; stp = 300;
        elseif sum(ismember(fish, 7)) > 0
          strt = 0; stp = 0;
        else
          strt = 0; stp = 1200;
        end
    end
    
    if id(2) == 14
        if sum(ismember(fish, 7)) > 0
          strt = 300; stp = 1000;
        else
          strt = 0; stp = 1000;
        end
    end    
    

end







end

