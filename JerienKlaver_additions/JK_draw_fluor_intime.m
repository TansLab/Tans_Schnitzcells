
function JK_draw_data_intime(change_over_time_all,branchdata)
color = ['r','g','b','c','m','y','k'];
count = 1;
figure;
for k = 1: length(branchdata)
    for i = branchdata(2,k)
        if (count>=8)
            count = 1
        end
            A = change_over_time_all{1,i};
            Ax = A(:,1);
            Ay = A(:,2);
            plot (Ax,Ay,color(count))
            xlabel('Time (mins)')
            ylabel('data')
            title('data of branch during time')
            hold on
            count = count+1
    end
end 
end
            