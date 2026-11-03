function [P,N,E,vc]=r30b_kinematics_many(s,z,h,t,p,A,tooth,side)
% 同一时刻逐空间积分点的坐标、目标法向、局部滑动方向。
s=s(:);z=z(:);h=h(:);[x,y,nx,ny]=r30b_surface(s,tooth,side,p,A);
th=A.phi0+A.omegaG*t;ct=cos(th);st=sin(th);
P=[(x+h.*nx)*ct-(y+h.*ny)*st,(x+h.*nx)*st+(y+h.*ny)*ct,z-(A.feedStart+A.feed*t)];
N=[nx*ct-ny*st,nx*st+ny*ct,zeros(size(s))];
axis=repmat(A.axisW(:)',numel(s),1);
vr=cross(A.omegaW*axis,P-repmat([A.center 0 0],numel(s),1),2)-...
 [-A.omegaG*P(:,2),A.omegaG*P(:,1),zeros(size(s))]+repmat([0 0 A.feed],numel(s),1);
E=vr-sum(vr.*N,2).*N;vc=sqrt(sum(E.^2,2));E=E./vc;
end
