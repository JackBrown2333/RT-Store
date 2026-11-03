function possible=r30b_gate(P,face,ta,tb,normalDepth)
% Conservative enclosing-cylinder rejection over the WHOLE time interval.
% Bounds every target-flank point plus h in [-normalDepth,stock]. No old
% contact masks, coarse point samples or assumed engagement phase are used.
p=P.p;A=P.installation;stock=P.cfg.stock_um/1000;
rb=p.calc.gear.baseDiameter_mm/2;ra=p.calc.gear.tipDiameter_mm/2;rm=(rb+ra)/2;
pressure=p.gear.pressure_deg*pi/180;
H0=p.calc.gear.pitchToothThickness_mm/(2*p.calc.gear.pitchRadius_mm)+tan(pressure)-pressure;
u=sqrt((ra/rb)^2-1);H1=H0-(u-atan(u));g=2*face(2)-3;
psi=(face(1)-1)*A.toothPitch+g*(H0+H1)/2;
rxy=(ra-rb)/2+ra*abs(H0-H1)/2+max(normalDepth,stock);
halfz=p.gear.faceWidth_mm/2;tm=(ta+tb)/2;radiusMotion=rm*abs(A.omegaG)*(tb-ta)/2;
zMotion=abs(A.feed)*(tb-ta)/2;angle=psi+A.phi0+A.omegaG*tm;
x=rm*cos(angle);y=rm*sin(angle);z=-(A.feedStart+A.feed*tm);
dx=max(0,abs(A.center-x)-rxy-radiusMotion);
dy=max(0,abs(A.axisY(2)*y+A.axisY(3)*z)-...
 abs(A.axisY(2))*(rxy+radiusMotion)-abs(A.axisY(3))*(halfz+zMotion));
axialLower=max(0,abs(A.axisW(2)*y+A.axisW(3)*z)-...
 abs(A.axisW(2))*(rxy+radiusMotion)-abs(A.axisW(3))*(halfz+zMotion));
possible=hypot(dx,dy)<=A.tipNew+1e-10&&axialLower<=A.axialLength/2+1e-10;
end
