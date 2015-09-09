function view_tracker(I, bbs_raw, bbs, indvl_detectors, frmInd)
figure(1); 
clf;
M = length(indvl_detectors); 

%define colours
col = uniqueColors(1, 6); 
col = col([1 2 4:6], :);

subplot(2,1,1); 
im(I);
colorbar off; 
title(num2str(frmInd));
draw_cat_id(bbs_raw(:, [1:5 6]), col, [], [], [], bbs_raw(:,6));

subplot(2,1,2); 
im(I);
colorbar off; 
draw_cat_id(bbs(bbs(:, 6) ~= 1, [1:4 6]), col, [], [], [], bbs(bbs(:, 6) ~= 1, 6));

end


function hs = draw_cat_id( bb, col, lw, ls, prop, ids )
% Draw single or multiple bbs to image (calls rectangle()).
%
% To draw bbs aligned with pixel boundaries, subtract .5 from the x and y
% coordinates (since pixel centers are located at integer locations).
%
% USAGE
%  hs = bbApply( 'draw', bb, [col], [lw], [ls], [prop], [ids] )
%
% INPUTS
%  bb     - [nx4] standard bbs or [nx5] weighted bbs
%  col    - ['g'] color or [kx1] array of colors
%  lw     - [2] LineWidth for rectangle
%  ls     - ['-'] LineStyle for rectangle
%  prop   - [] other properties for rectangle
%  ids    - [ones(1,n)] id in [1,k] for each bb into colors array
%
% OUTPUT
%  hs     - [nx1] handles to drawn rectangles (and labels)
%
% EXAMPLE
%  im(rand(3)); bbApply('draw',[1.5 1.5 1 1 .5],'g');
%
% See also bbApply, bbApply>embed, rectangle
[n,m]=size(bb); if(n==0), hs=[]; return; end
if(nargin<2 || isempty(col)), col=[]; end
if(nargin<3 || isempty(lw)), lw=2; end
if(nargin<4 || isempty(ls)), ls='-'; end
if(nargin<5 || isempty(prop)), prop={}; end
if(nargin<6 || isempty(ids)), ids=ones(1,n);  pos=0; end
% prepare display properties
prop{1}=['LineStyle' ls {} 'EdgeColor'];
prop{2}=['LineStyle' ls {} 'EdgeColor'];
tProp={'FontSize',10,'color','w','FontWeight','bold',...
  'VerticalAlignment','bottom'}; k=max(ids);
if(isempty(col)), if(k==1), col='g'; else col=hsv(k); end; end
hs=zeros(1,n);
% draw rectangles and optionally labels
for b=1:n, 
    if(ids(b)==1), col1='g'; else col1=col; end; 
    if(bb(b,5)<100),lw=2; else if(ids(b)==1), lw=10; else lw=4; end; end; 
    hs(b)=rectangle('Position',bb(b,1:4),prop{~(ids(b)==1)+1}{:},col1(mod(ids(b)-2,size(col1,1))+1,:),'LineWidth',lw); 

end
if(m==4), return; end; hs=[hs zeros(1,n)];
for b=1:n, if(ids(b)~=1), hs(b+n)=text(bb(b,1),bb(b,2),num2str(bb(b,5),4),tProp{:}); end; end
if(m==5), return; end;
% tProp={'FontSize',24,'color','w','FontWeight','bold',...
%   'VerticalAlignment','bottom'}; pos=1;
for b=1:n, if(ids(b)==1), hs(b+n)=text(bb(b,1),bb(b,2)+(bb(b,4)+30),num2str(bb(b,5),4),tProp{:}); end; end
for b=1:n, if(ids(b)~=1), hs(b+n)=text(bb(b,1)+bb(b,3)-2,bb(b,2),num2str(bb(b,6),4),tProp{:}); end; end
end