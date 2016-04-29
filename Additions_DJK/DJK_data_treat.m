%   [P D E G]= DJK_data_treat(dirname,extention);
%
% copied from data_treat, such that will work with DJK_tracker_original
%
% calculates the P,D,E,G frame-label-->next_frame-label correspondencise
% out of the raw output (.mch files) of the softassign-matching code.
% data will be read from dirname/*.extension (*.mch by default) files;
%
% P,D,E,G -are  correspondencies between label IDs 
% on consequtive frames matching same cells:
% P{nframe}(j) corresponds to the label from (nframe+1) frame 
% matching  the cell labeled j on the nframes-th frame, if it did not
% split.
% D,E -same for the case where there was a split into two other cells,
% G-logical index array indicating labels of 'ghost' cells which don't have
% a match in previous frames.
% T.B. 04/05

% JCR mods to input arguments & code to find match files using track range

function [P, D, E, G]= DJK_data_treat(p);

%%
P={};
D={};
E={};

i=1;
%% Loop over all frames 
for frameNum = p.manualRange(2:end)
    
    %% Administration (what is previous frame)
    mynum_t= str3(frameNum);
    yesterdayFrameNum = p.manualRange(find(p.manualRange==frameNum)-1);
    mynum_y= str3(yesterdayFrameNum);

    %% Load tracking file for frames n and n-1 into variable
    trackOutputFile = [p.tracksDir,p.movieName,'-djk-output-',...
                       mynum_y,'-to-',mynum_t,'.txt'];                   
    sizeCurrentP = load( trackOutputFile ,'ASCII');

    %% MW - unclear why this is needed, is data ever not Nx4 matrix?
    sz=size(sizeCurrentP(:),1)/4;      % sz = # all elements divided by four
    sizeCurrentP=reshape(sizeCurrentP,[sz,4]);    % reshapes matrix that was already Nx4 to Nx4?

    %% MW - after manual edit, data can be sorted in wrong order, fix that
    % here
    msorted = sortrows(sizeCurrentP,4);
    if ismember(0,msorted == sizeCurrentP)
        disp(['WARNING, last column in ' trackOutputFile ' was not sorted correctly. I have corrected this.']);
        sizeCurrentP = msorted;
    end
    
    %% MW TODO - make more clear what these labels stand for
    P1=sizeCurrentP(:,1); % MW nomencl. pn
    P2=sizeCurrentP(:,4); % MW nomencl. ch
    D1=sizeCurrentP(:,2); % MW nomencl. p1
    D2=sizeCurrentP(:,3); % MW nomencl. p2
   
    %P2sort=sort( round(P2),'ascend' ) ;

    %% EDIT MW (TODO: remove these comments if proven to work smoothly)
    % Data is now sorted if it was not to begin with (see above)
    % Old code:
    % First statement checks whether P2 is sorted 
    % Second checks whether first value is a one    
    % if( sum(P2==P2sort)~=size(P2,1) | P2sort(1)<1)
    % New code:
    if( P2(1)<1)
        %error(sprintf('last column in:  ''%s'' is not 1,2,3...',...
        error(['First value last column should be a one, not the case in ', ...
                       trackOutputFile] );
    else
        % Check whether last column is 1,2,3, correct otherwise
        if( sum(P2==[1:size(P2,1)]')~=size(P2,1) );
            warning(['last column in ',trackOutputFile,...
                    ' is not 1,2,3,... adjusting']);
            
            % correcting for missing numbers in P2
            P2bis=[1:max(P2)]'; P1bis=zeros(size(P2bis)); 
            D1bis=zeros(size(P2bis)); D2bis=zeros(size(P2bis));
            ii=0;         
            while( ii+size(P2,1) <= size(P2bis,1) )
                l=(  P2  == P2bis( [ 1:size(P2,1) ] + ii )  ) ;
                lbig=[logical( zeros(ii,1) ); l];
     
                P1bis( lbig )=  P1(l);
                D1bis( lbig )=  D1(l);
                D2bis( lbig )=  D2(l);     
                % [P1bis D1bis D2bis P2bis]
          
                ii=ii+1;
            end           
     
            P1=P1bis;
            P2=P2bis;
            D1=D1bis;
            D2=D2bis;
           
        end
    end

    %%
    
    % Create vars
    p1=zeros(size(P2));
    d1=zeros(size(P2));
    d2=zeros(size(P2));
    
    % Note P{nframe}(j) corresponds to the label from (nframe+1) frame 
    % matching  the cell labeled j on the nframes-th frame, if it did not
    % split.
    p1(P1(P1>0)) = ...
        P2(P1>0); % P1>0 are cells that did not divide, P2 are the labels in the next frame
    d1(D1(D1>0)) = ...
        P2(D1>0); % D1>0 are cells that did divide, P2 are the labels in the next frame
    d2(D2(D2>0)) = ...
        P2(D2>0); % D2>0 are cells that did divide, P2 are the labels in the next frame

    P{i}=p1;
    D{i}=d1;
    E{i}=d2;
    
    % For the first frame, ghost cells are needed
    if(i==1)
        G{1}=( ones(size(P2))==1 );
    end
    % And also create ghost links for unlinked cells
    G{i+1}= ~(P1 > 0 | D1 > 0 | D2 > 0);

    maxPnum=max( [max(P1) max(D1) max(D2)] );
 
    %cut extra zeros
    sizeCurrentP=size(P{i},1);
       
    % for frames other than the starting frame
    if(i>1)
        %[sizeP maxPnum m i];
        
        if(maxPnum>sizePreviousP)
           
            error(sprintf(['The match output files for frames %s %s might ',...
                           'not be consecutive matches: inconsistent IDs'],...
                           mynum_y,mynum_t) );
        end
        
        P{i}(sizePreviousP+1)=0;
        D{i}(sizePreviousP+1)=0; 
        E{i}(sizePreviousP+1)=0; 
        
        P{i}(sizePreviousP+1:max([sizeCurrentP sizePreviousP+1]) )=[];
        D{i}(sizePreviousP+1:max([sizeCurrentP,sizePreviousP+1]) )=[];
        E{i}(sizePreviousP+1:max([sizeCurrentP,sizePreviousP+1]) )=[];          
    end
        
    sizePreviousP=sizeCurrentP;
    i = i+1;
end
