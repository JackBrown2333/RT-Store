function G=r30d_wheel_geometry(P,t,A)
% Same IDEAL implicit solid as r30b_level, in world coordinates.
% Analytic outward gradient of the active max branch. Not a measured profile.
% V3: explicit N-by-3 arithmetic and guarded branches, including N=0 and N=1.
% branch: 1=helical flank, 2=outer cylinder, 3=axial end, 4=root core.
assert(ismatrix(P)&&size(P,2)==3&&all(isfinite(P(:))));n=size(P,1);
if isscalar(t),t=repmat(t,n,1);else,t=t(:);end;assert(numel(t)==n);
assert(all(isfinite(t)));
if n==0
 G=struct('level',zeros(0,1),'branch',zeros(0,1),'branchGap',zeros(0,1),...
  'gradient',zeros(0,3),'outwardNormal',zeros(0,3),...
  'partialTimeDerivative',zeros(0,1),'rigidTransportResidual',zeros(0,1));
 return;
end
aw=A.axisW(:)';ay=A.axisY(:)';qx=A.center-P(:,1);qy=P*ay';ax=P*aw';r=hypot(qx,qy);
assert(all(r>1e-12),'R30D:Axis','Gradient undefined at wheel axis.');
phi=atan2(qy,qx);phase=mod(ax-A.helicalP*(phi-A.omegaW*t)+A.lead/2,A.lead)-A.lead/2;
rb=A.rPitch/hypot(1,A.rPitch*A.slope/A.helicalP);up=sqrt((A.rPitch/rb)^2-1);
u=sqrt((max(r,A.rRoot)/rb).^2-1);half=A.lead/4-A.helicalP*(u-atan(u)-up+atan(up));
values=[abs(phase)-half,r-A.tipNew,abs(ax)-A.axialLength/2];
[level,branch]=max(values,[],2);sorted=sort(values,2,'descend');gap=sorted(:,1)-sorted(:,2);
core=r<=A.rRoot;cv=[r-A.rRoot,abs(ax)-A.axialLength/2];[cl,cb]=max(cv,[],2);
level(core)=cl(core);branch(core)=4;branch(core&cb==2)=3;gap(core)=abs(cv(core,1)-cv(core,2));
gr=(qx*[-1 0 0]+qy*ay)./repmat(r,1,3);
gp=(qx*ay+[qy zeros(n,2)])./repmat(r.^2,1,3);
gradient=zeros(n,3);gt=zeros(n,1);flank=branch==1;
derivative=-A.helicalP*sqrt(max(0,r.^2-rb^2))./(rb*r);
if any(flank)
 sg=reshape(sign(phase(flank)),[],1);dd=reshape(derivative(flank),[],1);
 gradient(flank,:)=sg*aw-A.helicalP*repmat(sg,1,3).*gp(flank,:)-repmat(dd,1,3).*gr(flank,:);
 gt(flank)=sg*A.helicalP*A.omegaW;
end
radial=branch==2|branch==4;gradient(radial,:)=gr(radial,:);
ends=branch==3;
if any(ends),gradient(ends,:)=reshape(sign(ax(ends)),[],1)*aw;end
gn=sqrt(sum(gradient.^2,2));assert(all(gn>0));normal=gradient./repmat(gn,1,3);
vw=cross(repmat(A.omegaW*aw,n,1),P-repmat([A.center 0 0],n,1),2);
G=struct('level',level,'branch',branch,'branchGap',gap,'gradient',gradient,...
 'outwardNormal',normal,'partialTimeDerivative',gt,'rigidTransportResidual',gt+sum(gradient.*vw,2));
end
