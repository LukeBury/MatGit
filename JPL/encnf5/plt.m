function h=plt(rr,proj,c,c3D)
% function h=plt(rr,proj,c,c3D)
% USE: plt0.m, original plt.m function.
% IN : rr = matrix of 2D, 3D row vectors to plot.
%      proj:1  xy
%           2  xz
%           3  yz
%           4  xyz
%	    5  all 4 projections.
%      c  = color & symbol, 'k','b+', etc.
%      c3D= color of 3D plt in proj 5, default='k'.

if nargin < 2
   proj=1;
   c='k';
elseif nargin < 3
   c='k';
elseif nargin < 4
   c3D='k';
end
if proj <= 1
   h=plot(rr(:,1),rr(:,2),c);
   xlabel('X');ylabel('Y');
elseif proj==2
   h=plot(rr(:,1),rr(:,3),c);
   xlabel('X');ylabel('Z');
elseif proj==3
   h=plot(rr(:,2),rr(:,3),c);
   xlabel('Y');ylabel('Z');
elseif proj==4
   h=plot3(rr(:,1),rr(:,2),rr(:,3),c);
   xlabel('X');ylabel('Y');zlabel('Z');axis('equal');
elseif proj==5
   h=subplot(2,2,1); plot(rr(:,1),rr(:,2),c);axis('equal');grid;
   xlabel('X');ylabel('Y');
   hold on;
   subplot(2,2,2); plot(rr(:,1),rr(:,3),c);axis('equal');grid;
   xlabel('X');ylabel('Z');
   hold on;
   subplot(2,2,3); plot(rr(:,2),rr(:,3),c);axis('equal');grid;
   xlabel('Y');ylabel('Z');
   hold on;
   subplot(2,2,4); plot3(rr(:,1),rr(:,2),rr(:,3),c3D);
   xlabel('X');ylabel('Y');zlabel('Z');axis('equal');grid;
   hold on;
end
axis('equal');
