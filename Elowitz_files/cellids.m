% [cellid PP DD EE]=      cellids(P,D,E,G)
%  used in recalc_schnitz()
% calculates the label_id->cellnumber correspondence
% 'cellid' array ( indexed: <cell id> = cellid{frame}(label) )
% and parent (PP), daughter one (DD), and daughter two (EE)arrays
% labeled: <cell id>= <??>(cell_id)
% T.B. 04/05 

function   [cellid PP DD EE]=      cellids(P,D,E,G);

cellid{1}=[1:size(P{1},1)]';

PP(1:size(P{1},1) )=0;
DD(1:size(D{1},1) )=0;%ST
EE(1:size(E{1},1) )=0;%ST
nextid=size(P{1},1);

for i=2:size(P,2)+1;
    k=1:size(P{i-1},1);
    
    
    l=P{i-1}>0;
    
    cellid{i}(P{i-1}(P{i-1}>0))=cellid{i-1}(P{i-1}>0);
    
    if(max(D{i-1}>0)>0)
    cellid{i}( D{i-1}( D{i-1}>0 ) ) =[nextid+1:nextid+ sum(D{i-1}>0) ];
    PP([nextid+1:nextid+ sum( D{i-1}>0 ) ])= cellid{i-1}( D{i-1}>0 ) ;
    nextid=nextid+ sum(D{i-1}>0);
    end
    
    if(max(E{i-1}>0)>0)
      
    cellid{i}( E{i-1}( E{i-1}>0 ) ) =[nextid+1:nextid+ sum(E{i-1}>0)];
    PP([nextid+1:nextid+ sum(E{i-1}>0)])= cellid{i-1}( E{i-1}>0 );
    nextid=nextid+ sum(E{i-1}>0);
    end
    
  
    
         if(sum(G{i})>0) 
   cellid{i}( G{i}>0 ) =[nextid+1:nextid+ sum(G{i}>0)];
    PP(  [nextid+1 : nextid + sum( G{i} ) ] )=0;
    nextid=nextid+ sum(G{i}>0);
   end

    
   if( max(D{i-1}>0) >0)
       
    lD=(D{i-1}>0 );
    
 %   lD=(cellid{i-1}(lD)>0);

%    cellid{i-1}(lD)
%    if(i==35)
%       cellid{i}(D{i-1}(lD));
%       G{i-1}(lD)'
%       cellid{i-1}(lD)
%       
%       DD(cellid{i-1}(lD)) =  cellid{i}( D{i-1}( lD ) );
%    end
     DD(cellid{i-1}(lD))=cellid{i}( D{i-1}( lD ) );

   
   end
   if( max(E{i-1}>0) >0)
     EE(cellid{i-1}( E{i-1}>0))...
                                    =cellid{i}( E{i-1}( E{i-1}>0 ) );
   end
end

%keyboard;

if( max(size(DD)) <nextid)
DD(nextid)=0;
end

if( max(size(EE)) <nextid)
EE(nextid)=0;
end

if( max(size(PP)) <nextid)
PP(nextid)=0;
end















