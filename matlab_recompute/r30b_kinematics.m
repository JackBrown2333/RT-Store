function K=r30b_kinematics(p,A,t,tooth,side,s,z,height)
% 当前实际三维接触点、法线、滑动速度，全部来自同一t的砂轮/工件运动。
[x,y,nx,ny,tx,ty]=r30b_surface(s,tooth,side,p,A);
th=A.phi0+A.omegaG*t;ct=cos(th);st=sin(th);
P=[(x+height*nx)*ct-(y+height*ny)*st;...
   (x+height*nx)*st+(y+height*ny)*ct;z-(A.feedStart+A.feed*t)];
n=[nx*ct-ny*st;nx*st+ny*ct;0];
es=[tx*ct-ty*st;tx*st+ty*ct;0];ez=[0;0;1];
vTool=cross(A.omegaW*A.axisW(:),P-[A.center;0;0]);
vWork=cross([0;0;A.omegaG],P-[0;0;P(3)])-[0;0;A.feed];
vr=vTool-vWork;vTang=vr-dot(vr,n)*n;
vc=norm(vTang);
if vc<1e-8,error('R30:Velocity','局部滑动速度退化，不能计算磨削力。');end
slide=vTang/vc;
K=struct('point',P,'normal',n,'slide',slide,'es',dot(slide,es),...
    'ez',dot(slide,ez),'vc',vc,'vf',abs(A.feed));
end
