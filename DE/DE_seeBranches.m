function br = DE_seeBranches(p, schnitzcells, fitTime)

close all;
N = 3;
%fitTime = [3 2000];

%'all'
%'fitTime'
%'fitTime_cycle'

what_filter = 'all';

% filter schnitzcells by time-related properties:
switch what_filter
    case 'all'
        disp('filter : none')
        s = DJK_selSchitzesToPlot(schnitzcells, 'P', @(x) 1);    
    case 'fitTime'
        disp('filter : fitTime')
        % restrict to fitTime
        s = DJK_selSchitzesToPlot(schnitzcells, 'time', @(x) x(1) > fitTime(1) & x(1) < fitTime(2)); 
    case 'cycle'
        % restrict to completeCycle
        s = DJK_selSchitzesToPlot(schnitzcells, 'completeCycle', @(x) x ~= 0);
    case 'fitTime_cycle'
        disp('filter : fitTime_cycle')
        s = DJK_selSchitzesToPlot(schnitzcells, 'time', @(x) x(1) > fitTime(1) & x(1) < fitTime(2)); 
        s = DJK_selSchitzesToPlot(s, 'completeCycle', @(x) x ~= 0);
end


% if no R time vector is present (ony G-time is there), create R-time
% same as G:
if ~isfield(s, 'dR5_time')
    for i=1:length(s)
        s(i).dR5_time =  s(i).dG5_time;
        s(i).R_time =  s(i).G_time;
        s(i).dR5 = s(i).dG5;
        s(i).dR5_cycCor = s(i).dG5_cycCor;
        s(i).R5_mean = s(i).G5_mean;
        s(i).R5_mean_cycCor = s(i).G5_mean_cycCor;
        s(i).muP15_fitNew_atdR5 = s(i).muP15_fitNew_atdG5;
        s(i).muP15_fitNew_atdR5_cycCor = s(i).muP15_fitNew_atdG5_cycCor;
    end

end

br = DJK_getBranches(p,s,'dataFields',...
    {'dR5_time'  'R_time' ...
    'dR5' 'dR5_cycCor' 'R5_mean' 'R5_mean_cycCor' ...
    'muP15_fitNew_atdR5' 'muP15_fitNew_atdR5_cycCor', 'muP15_fitNew'...
    'frames'}, ...
    'fitTime', fitTime);






fx1 = 'R_time';
fy1 = 'R5_mean';

fx2 = 'R_time';
fy2 = 'muP15_fitNew';


for i=1:length(br)
    % number of all schnitzcells belonging to this branch like so: [1 1 1 2 2 2 2 2 4 4 4 4 7 7 7 7]
     nr = br(i).schnitzNrs;
     % find the differencial:
     c = diff(nr);
     % compensate for one lost entry:
     c = [c(1),c];
     c = double(logical(c));
     
     % this matrix has 1 when a new cell was born:
     br(i).division_events=[br(i).schnitzNrs; br(i).R_time; c;  br(i).frames];
end



h1=figure; hold on;
h2=figure; hold on;

for i=1:length(br)
    figure(h1);
    plot(br(i).(fx1),br(i).(fy1),'-','color',[0.5 0.5 0.5])
end

colors = jet(N);
selected_branches = [];

for i=1:N
    
    figure(h1);
    current_color = colors(i,:);
   % w = waitforbuttonpress;
    %key1 = double(get(gcf, 'CurrentCharacter'));
    
    [x,y] = ginput(1);
    clicked_branch = [];
    min_distance = 1000;

    for k = 1:length(br)
        dx = br(k).(fx1) - x;
        dy = br(k).(fy1) - y;
        ds = sqrt(dx.^2+dy.^2);
        br(k).ds_min = min(ds);
    end
    
    all_min_distances = [br.ds_min]';
    clicked_branch = find(all_min_distances == min(all_min_distances))
    
    for m = 1:length(clicked_branch)
        
        branch_m = br(clicked_branch(m));
        
        clicked_branch(m);
        divisions_m = logical(branch_m.division_events(3,:));

        
        figure(h1);
        plot(branch_m.(fx1),branch_m.(fy1),'--','color',current_color,'LineWidth',2)
        
        plot(branch_m.(fx1)(divisions_m),branch_m.(fy1)(divisions_m),'rx','MarkerSize',8)%,'color',current_color,'LineWidth',2)

        title([fy1 ' vs time. Branch number ' num2str(clicked_branch(m)) ])
        
        
        figure(h2);
        plot(branch_m.(fx2),branch_m.(fy2),'--','color',current_color,'LineWidth',2)
        title([fy2  ' vs time.'])
        
        selected_branches = [selected_branches; i*ones(size(clicked_branch)) clicked_branch ];
    end   
end



selected_branches

























if 0


[~,ind]=unique(taken,'rows');
large_fluct=large_fluct(ind);

bignum=100;
colors=jet(bignum);


if tr>0
    
    N_temp=floor(bignum/length(large_fluct));

    colors_temp=colors(N_temp*[1:1:length(large_fluct)],:);

    
    HH1=[];
    HHc1=[];
    HHx1=[];
    figure(h1);
    
    
    for j=1:length(large_fluct)
        b_j=br(large_fluct(j));
        %all_colors=[all_colors;rand(1,3)];
        HH1(j)=plot(b_j.(fx1),b_j.(fy1),'-','color', colors_temp(j,:),'LineWidth',2);
        
        %plot frame markers
        HHx1(j)=plot(b_j.(fx1),b_j.(fy1),'kx','MarkerSize',6,'LineWidth',2);        
        
        % index of the "new cell"- event
        [~,new_cell_ind]=find(b_j.new_cell_events(3,:)==1);
        HHc1(j)=plot(b_j.(fx1)(new_cell_ind),b_j.(fy1)(new_cell_ind),...
            'ko','MarkerSize',10,'MarkerFaceColor', colors_temp(j,:),'MarkerEdgeColor','r','LineWidth',2);
        
              
        for kk=1:length(b_j.new_cell_events)
            if b_j.new_cell_events(3,kk)==1
                % name of the schnitz  + its frame:
                name_and_frame=num2str([b_j.new_cell_events(1,kk), b_j.new_cell_events(4,kk)]);
                text(b_j.(fx1)(kk),b_j.(fy1)(kk),name_and_frame,'FontSize',15)
            end
        end
        
    end
    
    legend(HH1,num2str(large_fluct));
    
    fx1=strrep(fx1,'_','-');
    fy1=strrep(fy1,'_','-');
    title([pos_name '; ' fy1 ' vs ' fx1])
    
    
    
    
    
    %now plot the same branches for the other fields:
    HH2=[];
    HHc2=[];
    HHx2=[];
    figure(h2);
    
    for j=1:length(large_fluct)
        
        b_j=br(large_fluct(j));

        HH2(j)=plot(br(large_fluct(j)).(fx2),br(large_fluct(j)).(fy2),'-','color', colors_temp(j,:),'LineWidth',2);
        
        
        %plot frame markers
        HHx2(j)=plot(b_j.(fx2),b_j.(fy2),'kx','MarkerSize',6,'LineWidth',2);        
        
        % index of the "new cell"- event
        [~,new_cell_ind]=find(b_j.new_cell_events(3,:)==1);
        HHc2(j)=plot(b_j.(fx2)(new_cell_ind),b_j.(fy2)(new_cell_ind),...
            'ko','MarkerSize',10,'MarkerFaceColor', colors_temp(j,:),'MarkerEdgeColor','r','LineWidth',2);
        
        j
        for kk=1:length(b_j.new_cell_events)
            if b_j.new_cell_events(3,kk)==1
                % name of the schnitz  + its frame:
                name_and_frame=num2str([b_j.new_cell_events(1,kk), b_j.new_cell_events(4,kk)]);
                text(b_j.(fx2)(kk),b_j.(fy2)(kk),name_and_frame,'FontSize',15)
            end
        end
        
    end
        
    fx2=strrep(fx2,'_','-');
    fy2=strrep(fy2,'_','-');
    title([pos_name '; ' fy2 ' vs ' fx2])
    
else %if Tr=0 - no colors needed
    figure(h1);
    fx1=strrep(fx1,'_','-');
    fy1=strrep(fy1,'_','-');
    title([pos_name '; ' fy1 ' vs ' fx1])
    
    figure(h2);
    fx2=strrep(fx2,'_','-');
    fy2=strrep(fy2,'_','-');
    title([pos_name '; ' fy2 ' vs ' fx2])
end
%legend(HHc);

end