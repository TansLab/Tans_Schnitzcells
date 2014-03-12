% to do regression for both parent and switching schnitz then find the crossing

regX=[[AllData{1,2}{1,1} AllData{1,2}{1,end}]]
regY=[AllData{1,3}{1,1} AllData{1,3}{1,end}]

plot(regX,regY,'-')


done = 0;
while ~done
    
    datacursormode on;
    
    %first regression :
    waitforbuttonpress;
    dcm_obj = datacursormode(gcf);
    c_info = getCursorInfo(dcm_obj);
    Index1 = c_info.DataIndex;
    waitforbuttonpress;
    dcm_obj = datacursormode(gcf);
    c_info = getCursorInfo(dcm_obj);
    Index2 = c_info.DataIndex;
   
    regX([Index1:Index2]);
    regY([Index1:Index2]);
    fitpar1 = polyfit(regX([Index1:Index2]),regY([Index1:Index2]),1) 
    fitpar2 = fitpar * 10^3
    fitfunct = polyval(fitpar1,regX([1:Index2]));
    hold on;
    plot(regX([1:Index2]),fitfunct);
    hold on;
    reply1 = input('correct (y/n)', 's');
    reply1
   
    %second regression:
    waitforbuttonpress;
    dcm_obj = datacursormode(gcf);
    c_info = getCursorInfo(dcm_obj);
    Index3 = c_info.DataIndex;
    waitforbuttonpress;
    dcm_obj = datacursormode(gcf);
    c_info = getCursorInfo(dcm_obj);
    Index4 = c_info.DataIndex;
    
    regX([Index3:Index4]);
    regY([Index3:Index4]);
    fitpar3 = polyfit(regX([Index3:Index4]),regY([Index3:Index4]),1) 
    fitpar4 = fitpar * 10^3
    fitfunct = polyval(fitpar3,regX([1:Index4]));
    hold on;
    plot(regX([1:Index4]),fitfunct);
    hold on;
    reply2 = input('correct (y/n)', 's');
    reply2
    
    if and(reply1=='y',reply2=='y')
        done = 1;
    end
   
    Tsw=(fitpar3(2)-fitpar1(2))/(fitpar1(1)-fitpar3(1))
    Tsw
    
    
end
hold off;
close(gcf);