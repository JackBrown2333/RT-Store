function [x,y,nx,ny,tx,ty]=r30b_surface(s,tooth,side,p,A)
% 基圆以上渐开线；s为齿廓弧长，单位mm。齿根圆角不包含在内。
tooth=double(tooth);side=double(side);s=double(s);
rb=p.calc.gear.baseDiameter_mm/2;r=sqrt(rb^2+2*rb*s);
u=sqrt(2*s/rb);alpha=atan(u);a=p.gear.pressure_deg*pi/180;
H=p.calc.gear.pitchToothThickness_mm/(2*p.calc.gear.pitchRadius_mm)+tan(a)-a-(u-alpha);
g=2*side-3;psi=(tooth-1)*A.toothPitch+g*H;heading=psi-g*alpha;
x=r.*cos(psi);y=r.*sin(psi);tx=cos(heading);ty=sin(heading);
nx=-g*ty;ny=g*tx;
end
