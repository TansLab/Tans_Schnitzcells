function [oneout,twoout]=padandshift(one,two,trans,len,samesize);
% function [onepad,twopadshift]=padandshift(one,two,trans,len,samesize);
%
% returns image one padded, and image two padded and shifted.

if nargin<4
    len=0;
    samesize=0;
elseif nargin<5
    samesize=0;
end

padsize= floor(max(abs(trans))+len)+1;
% today
onepad= zeros(size(one) + 2*padsize);
onepad(padsize + 1:end - padsize, padsize + 1: end - padsize)= one;
% yesterday
twopadshift= zeros(size(two) + 2*padsize);
twopadshift(padsize + 1 - trans(1):end - padsize - trans(1),...
    padsize + 1 - trans(2): end - padsize - trans(2))= two;

if samesize
    sizei=max(size(onepad,1),size(twopadshift,1));
    sizej=max(size(onepad,2),size(twopadshift,2));
    oneout=zeros(sizei,sizej);
    oneout(1:size(onepad,1),1:size(onepad,2))=onepad;
    twoout=zeros(sizei,sizej);
    twoout(1:size(twopadshift,1),1:size(twopadshift,2))=twopadshift;
else
    oneout=onepad;
    twoout=twopadshift;
end