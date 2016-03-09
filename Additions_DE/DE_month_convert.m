function month_out=DE_month_convert(month_in)

% Routine that converts 'Month' into number and number into 'Month'.
% Used in createMetaData. Related to time extracting from entries of kind
% '22-1-2014'...


if ~ischar(month_in)
    
    switch month_in
        case 1
            month_out='Jan';
        case 2
            month_out='Feb';
        case 3
            month_out='Mar';
        case 4
            month_out='Apr';
        case 5
            month_out='May';
        case 6
            month_out='Jun';
        case 7
            month_out='Jul';
        case 8
            month_out='Aug';
        case 9
            month_out='Sep';
        case 10
            month_out='Oct';
        case 11
            month_out='Nov';
        case 12
            month_out='Dec';               
    end
    
    
elseif ischar(month_in)
    
     switch month_in
        case 'Jan'
            month_out=1;
        case 'Feb'
            month_out=2;
        case 'Mar'
            month_out=3;
        case 'Apr'
            month_out=4;
        case 'May'
            month_out=5;
        case 'Jun'
            month_out=6;
        case 'Jul'
            month_out=7;
        case 'Aug'
            month_out=8;
        case 'Sep'
            month_out=9;
        case 'Oct'
            month_out=10;
        case 'Nov'
            month_out=11;
        case 'Dec'
            month_out=12;               
    end
          
end
