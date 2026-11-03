function Q=r30b_init(P,o,face)
qopts=struct('ns',o.ns,'nz',o.nz,'order',o.order,'z0',-P.p.gear.faceWidth_mm/2,...
 'z1',P.p.gear.faceWidth_mm/2,'tooth',face(1),'side',face(2));
Q=r30b_quadrature(P,qopts);Q.passStart=Q.h;Q.inPass=false;Q.passId=0;
sc=(Q.se(1:end-1)'+Q.se(2:end)')/2;
[x,y,nx,ny]=r30b_surface(sc,face(1),face(2),Q.p,Q.A);
Q.centerProfile=[x y nx ny];Q.zc=(Q.ze(1:end-1)+Q.ze(2:end))/2;
Q.rollStep=diff(sqrt(2*Q.se/Q.rb))';Q.historyVolume_mm3=0;
Q.weightsFine=diff(Q.support_s)*diff(Q.support_z)';
end
