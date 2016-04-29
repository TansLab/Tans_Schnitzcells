A=imread('F:\Datasets\test_rotate_72x5degrees\pos2\pos2-p-1-000.tif');
for i=1:71
    B=imrotate(A,5*i,'bilinear','crop');
    if i<10
        imwrite(B,['F:\Datasets\test_rotate_72x5degrees\pos2\pos2-p-1-00' num2str(i) '.tif']);
    elseif i<100
        imwrite(B,['F:\Datasets\test_rotate_72x5degrees\pos2\pos2-p-1-0' num2str(i) '.tif']);
    else
        imwrite(B,['F:\Datasets\test_rotate_72x5degrees\pos2\pos2-p-1-' num2str(i) '.tif']);
    end
end
%%
C=imread('F:\Datasets\test_rotate_72x5degrees\pos1crop\images\pos1crop-p-1-000.tif');
for i=1:71
    D=imrotate(C,5*i,'bilinear','crop');
    if i<10
        imwrite(D,['F:\Datasets\test_rotate_72x5degrees\pos1crop\images\pos1crop-p-1-00' num2str(i) '.tif']);
    elseif i<100
        imwrite(D,['F:\Datasets\test_rotate_72x5degrees\pos1crop\images\pos1crop-p-1-0' num2str(i) '.tif']);
    else
        imwrite(D,['F:\Datasets\test_rotate_72x5degrees\pos1crop\images\pos1crop-p-1-' num2str(i) '.tif']);
    end
end