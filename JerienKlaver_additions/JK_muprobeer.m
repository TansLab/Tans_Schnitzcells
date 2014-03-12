%JK 08-02-08 This function adds an array to each schnitz wich contains 
% growth rate data. To compute growth rate for a certain time point, 
% 'length' length-measurements are used (10 before and 10 after timeframe
% for wich growthrate is computed)


function [schnitzcells] = JK_getmuarray (schnitzcells,slidingwindowlength);
    
for j = 1:length(schnitzcells);

    muarray = [];
    x = schnitzcells(1,(j)).mins;
    y = schnitzcells(1,(j)).lengthMicrons;
    fit = polyfit(x, log(y), 1);
    Growth = [fit(1)*60/log(2)];
    
    if length(x)< slidingwindowlength+1;
     for i=1:length(x);
        muarray(i)= Growth;
     end
    

    else length(x)> slidingwindowlength+1;
        for i=1:(slidingwindowlength/2);
        muarray(i)= Growth;
        end

        for i=(slidingwindowlength/2+1:length(schnitzcells(1,(j)).mins)-slidingwindowlength/2);
            %this means: for middel part of array
     
        x = schnitzcells(1,(j)).mins(i-slidingwindowlength/2:i+slidingwindowlength/2); 
        y = schnitzcells(1,(j)).lengthMicrons(i-slidingwindowlength/2:i+slidingwindowlength/2);
        fit = polyfit(x, log(y), 1);
        muarray(i)=fit(1)*60/log(2);
        end
        
        for i=1:(slidingwindowlength/2);
        muarray(i)= muarray(slidingwindowlength/2+1);
        end
        
        for i=length(schnitzcells(1,j).mins)-slidingwindowlength/2:length(schnitzcells(1,j).mins);
        muarray(i)= muarray(schnitzcells(1,j).mins)-slidingwindowlength/2-1);
        end
    end
       
schnitzcells(1,j).Growth = muarray;
end
schnitzcells = orderfields(schnitzcells,{'mins' ;'lengthMicrons' ;'Growth' ;'MY' ;'B' ;'D' ;'E' ;'FC' ;'FCs' ; 
'FCsframes' ; 'FCsmins' ; 'FCsvols' ;
'FY' ; 'FY_interps_Ydot' ; 'FYs' ; 'FYsframes' ; 'FYsmins' ; 'FYsvols' ;  'LCs' ; 'LYs' ; 'MC' ; 
'MCs' ;  'MYs' ; 'N' ; 'P' ; 'SC' ; 'SY' ; 'SZ' ; 'ang' ; 'approved' ; 'cellno' ; 'cenX' ; 'cenY' ; 
'cenx' ; 'ceny' ; 'dFCdt' ; 'dFCdvol' ; 'dFCmins' ; 'dFYdt' ; 'dFYdvol' ; 'dFYmins' ; 'datetime' ; 
'frames' ; 'len' ;'medC' ; 'medY' ;  'phase' ; 'phase_interps_Cdot' ; 
'phase_interps_Ydot' ; 'vol_interps_Ydot' ; 'volume' ; 'wid' ; 'width_Microns'});
end