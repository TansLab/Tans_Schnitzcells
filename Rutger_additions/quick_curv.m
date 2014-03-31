function [curv1, curv2] = quick_curv(x,y)
%% Returns a the curvatures at each point of a curve given as (x,y). For the
% best use y is a smooth function of x.

% Basic idea is found somewhere in Matlab news (online).

% The idea is to calculate curvature in each point of a cell over the whole
% length (this function). Then create a field:
% schnitzcells.curvature = curve;
% Then when the whole structure Schnitzcells is assembled, one can sort
% cells by their curvatures. 

%%
dx = gradient(x);
dy = gradient(y);
curv1 = gradient(atan2(dy,dx)) ./ hypot(dx,dy);


dx  = gradient(x);
ddx = gradient(dx);
dy  = gradient(y);
ddy = gradient(dy);
num   = dx .* ddy - ddx .* dy;
denom = dx .* dx + dy .* dy;
denom = sqrt(denom);
denom = denom.^3;
curv2 = num ./ denom;
curv2(denom < 0) = NaN;

% Demo with a circle and cosine. Just copypaste the inside of the if structure:
if 0
    R = 2;
    step = 1e-5;
    x1 = [-R:step:R];
    y1 = -(R^2-x1.^2).^0.5;
    
    y2 = (R^2-x1.^2).^0.5;
    x2 = x1 + 2*R + step; %shift to the right
    x = [x1, x2];
    y = [y1, y2];
    [curv1 curv2] = quick_curv(x,y);
    subplot(1,2,1),plot(x,y); hold on
    subplot(1,2,1),plot(x,curv1,'r')
    subplot(1,2,1),plot(x,curv2,'g--')
    title({['Two half-circles. R = ' num2str(R)], ['The peaks in curvatures are due to incontinuity']})
    legend('Circles','Curvatures')
    
    % cos
    x = [-3*pi:step:3*pi];
    y = cos(x);
    [curv1 curv2] = quick_curv(x,y);
    subplot(1,2,2),plot(x,y); hold on
    subplot(1,2,2),plot(x,curv1,'r')
    subplot(1,2,2),plot(x,curv2,'g--')
    title(['Cosine'])
    legend('Cosine','Curvatures')
 end