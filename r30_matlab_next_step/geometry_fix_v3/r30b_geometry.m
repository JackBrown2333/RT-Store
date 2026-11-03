function [area,s,z,volume]=r30b_geometry(Q)
% Stable area/centroid of clipped numerical supports. Same thin-layer volume.
% Translate EACH polygon to its first vertex before the shoelace sums.
% Global-coordinate cancellation can move a thin sliver's centroid outside it.
if isempty(Q),area=zeros(0,1);s=area;z=area;volume=area;return;end
n=size(Q,1);m=size(Q,2);assert(m>=3&&size(Q,3)>=6&&all(isfinite(Q(:))),...
 'R30B:PolygonInput','Nonfinite or incomplete support polygon.');
ox=Q(:,1,1);oy=Q(:,1,2);
x=bsxfun(@minus,Q(:,:,1),ox);y=bsxfun(@minus,Q(:,:,2),oy);
xn=x(:,[2:end 1]);yn=y(:,[2:end 1]);cross=x.*yn-xn.*y;
den=sum(cross,2);area=abs(den)/2;
s=ox+mean(x,2);z=oy+mean(y,2);ok=den~=0;
sx=sum((x+xn).*cross,2);zy=sum((y+yn).*cross,2);
s(ok)=ox(ok)+sx(ok)./(3*den(ok));z(ok)=oy(ok)+zy(ok)./(3*den(ok));
% A clipped rectangle is convex; its area center must stay in its bounds.
tol=128*eps(max(1,max(abs(reshape(Q(:,:,1:2),[],1)))));
assert(all(isfinite(s))&&all(isfinite(z))&&...
 all(s>=min(Q(:,:,1),[],2)-tol&s<=max(Q(:,:,1),[],2)+tol)&&...
 all(z>=min(Q(:,:,2),[],2)-tol&z<=max(Q(:,:,2),[],2)+tol),...
 'R30B:PolygonCentroid','Clipped support centroid is outside its polygon bounds.');
volume=zeros(n,1);
for j=2:m-1
 ar=abs(x(:,j).*y(:,j+1)-x(:,j+1).*y(:,j))/2;
 volume=volume+ar.*(Q(:,1,6)+Q(:,j,6)+Q(:,j+1,6))/3;
end
end
