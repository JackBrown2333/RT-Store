function [area,s,z,volume]=r30b_geometry(Q)
% Area, true area centroid, and affine normal-depth volume of padded polygons.
n=size(Q,1);m=size(Q,2);x=Q(:,:,1);y=Q(:,:,2);
xn=x(:,[2:end,1]);yn=y(:,[2:end,1]);cross=x.*yn-xn.*y;den=sum(cross,2);area=abs(den)/2;
s=mean(x,2);z=mean(y,2);ok=abs(den)>1e-25;
sx=sum((x+xn).*cross,2);zy=sum((y+yn).*cross,2);
s(ok)=sx(ok)./(3*den(ok));z(ok)=zy(ok)./(3*den(ok));volume=zeros(n,1);
for j=2:m-1
 ar=abs((Q(:,j,1)-Q(:,1,1)).*(Q(:,j+1,2)-Q(:,1,2))-(Q(:,j,2)-Q(:,1,2)).*(Q(:,j+1,1)-Q(:,1,1)))/2;
 volume=volume+ar.*(Q(:,1,6)+Q(:,j,6)+Q(:,j+1,6))/3;
end
end
