function imagescaled = NW_rescalecenterimage_affine(origimage,rescalefactor)
% rescales image by factor 'rescalefactor' (2dim vector) along the diagonals.
% Then centers it so that new central point is the same as the original center.
% Then image is cut to have original dimensions
% works only for scaling factors >=1!

%for debugging: show figures
showfig=0;

if length(rescalefactor)==1
    disp(['Only one rescale factor given instead of two. Will perform homogeneous stretching along both diagonals']);
    rescalefactor=[rescalefactor, rescalefactor];
end
% stretch factors
s1=rescalefactor(1);
s2=rescalefactor(2);


% input image
%origimage=checkerboard(10,2); debugging
%origimage=imread('D:\DummyExp\2012-fl-uo\TestData\pos5crop-c-035.tif');
if showfig==1
    h1=figure();
    clf
    colormap('gray');
    %imshow(origimage,'Border','loose')
    imagesc(origimage)
    hold on
    axis on
end

%max coordinates
dims=size(origimage); 
d1=dims(1); d2=dims(2);
%diagonal vectors from start to end
% v1=top left -> bottom right. v2=bottom left -> top right
v1=[d2-1,d1-1]; %(1=initial coordinate in upper left corner (not =0!))
v2=[d2-1,1-d1];
% better readability (like in script)
x1=v1(1);y1=v1(2);x2=v2(1);y2=v2(2); %y2<0. x2>0

% transformation matrix + elements;
c=x1*y2-x2*y1; % pre constant
t11=s1*x1*y2-s2*x2*y1;
t12=y1*y2*(s1-s2);
t21=x1*x2*(s2-s1);
t22=s2*x1*y2-s1*x2*y1;
T= 1/c*[t11  t12;
    t21  t22;
    0     0]; % no shift

t_aff = maketform('affine',T);
image_affine = imtransform(origimage,t_aff,'FillValues',0.3,'XYScale',1);

if showfig==1
    h2=figure();
    hold off
    clf
    colormap('gray');
    %imshow(image_affine,'Border','loose')
    imagesc(image_affine)
    hold on
    axis on
    %axis image
end

%rescenter image and cut such that it has the former dimensions
imagedim=size(origimage);   % size of rescaled origimage must be cut to this initial size
              
% find centre points of image (not yet rounded) and their
% shift after rescaling
centerpoints=0.5*(size(origimage)-[1,1])+[1,1]; % add [1,1] at the end because first pixel in matlab has coordinate (1,1) and not (0,0). For calculating the shift of the centre, it is however subtracted again anyway and so could also be constantly(!) ignored
newmaxcoordinates=max(s1*[x1,y1]+[1,1],s2*[x2,-y2]+[1,1]); %max should be for same argument in 1st and 2nd coordinate
xnewmax=newmaxcoordinates(1); ynewmax=newmaxcoordinates(2);
newcenterpoints=[1,1]+0.5*([ynewmax,xnewmax]-[1,1]); % [row, column]=[y,x]!

%newcenterpoints=s1*(centerpoints-[1,1])+[1,1];
centerpointsshift=round(newcenterpoints-centerpoints);
rowshift=centerpointsshift(1); columnshift=centerpointsshift(2); % only for better readability;
              

% fit in old frame
imagescaled=zeros(size(origimage));
size(imagescaled);
size(image_affine);
1+rowshift;
imagedim(1)+rowshift;
1+columnshift;
imagedim(2)+columnshift;
imagescaled(1:imagedim(1),1:imagedim(2))=image_affine(1+rowshift:imagedim(1)+rowshift,1+columnshift:imagedim(2)+columnshift);
                        % how to choose the shift: imaginge scale factor=1, then rowshift, columnshift=0 
                        % and first coordinate must be mapped on first coordinate

                        % for some reason, this leads to the fact that (when scaling by a
                        % factor 2) the new center is a 2x2 mattrix with oldcentercoordinates and
                        % -1 in each coordinate, e.g. (45,50) -> (44 45, 49 50)
                        % should not be a big problem because origimage
                        % still will be shifted


if showfig==1
    h3=figure();
    hold off
    clf
    colormap('gray');
    %imshow(imagescaled,'Border','loose')
    imagesc(imagescaled)
    hold on
    axis on
    %axis image
end



















%T = [1  0.1; 
%     1  1;
%     0  0];
%t_aff = maketform('affine',T);
%I_affine = imtransform(I,t_aff,'FillValues',0.3);


