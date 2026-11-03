function v=r30b_level(s,z,h,t,p,A,tooth,side)
% 当前砂轮实体的隐式值；隐式值不是法向距离。
if nargin<7,tooth=17;end;if nargin<8,side=1;end
[x,y,nx,ny]=r30b_surface(s,tooth,side,p,A);
th=A.phi0+A.omegaG*t;ct=cos(th);st=sin(th);
X=(x+h.*nx)*ct-(y+h.*ny)*st;Y=(x+h.*nx)*st+(y+h.*ny)*ct;
Z=z-(A.feedStart+A.feed*t);
ax=A.axisW(2)*Y+A.axisW(3)*Z;qx=A.center-X;qy=A.axisY(2)*Y+A.axisY(3)*Z;
r=hypot(qx,qy);eta=atan2(qy,qx)-A.omegaW*t;
phase=mod(ax-A.helicalP*eta+A.lead/2,A.lead)-A.lead/2;
half=r30b_worm_width(max(r,A.rRoot),A);
v=max(max(abs(phase)-half,r-A.tipNew),abs(ax)-A.axialLength/2);
core=r<=A.rRoot;v(core)=max(r(core)-A.rRoot,abs(ax(core))-A.axialLength/2);
end
