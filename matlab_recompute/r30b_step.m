function [Q,volume,nCandidates,change]=r30b_step(Q,t,o,screen)
% A fixed material column is updated once; the same positive weights give V.
if nargin<4,screen=true;end
H=reshape(Q.h,Q.ns*Q.nz,Q.n^2);A=Q.A;cp=Q.centerProfile;
change=struct('t_s',t,'index',zeros(0,1),'before_mm',zeros(0,1),'after_mm',zeros(0,1));
if screen
 lo=reshape(min(H,[],2),Q.ns,Q.nz);hi=reshape(max(H,[],2),Q.ns,Q.nz);hc=(lo+hi)/2;
 th=A.phi0+A.omegaG*t;
 X=(cp(:,1)+hc.*cp(:,3))*cos(th)-(cp(:,2)+hc.*cp(:,4))*sin(th);
 Y=(cp(:,1)+hc.*cp(:,3))*sin(th)+(cp(:,2)+hc.*cp(:,4))*cos(th);
 Z=Q.zc-(A.feedStart+A.feed*t);ax=A.axisW(2)*Y+A.axisW(3)*Z;
 qx=A.center-X;qy=A.axisY(2)*Y+A.axisY(3)*Z;r=hypot(qx,qy);
 phase=mod(ax-A.helicalP*(atan2(qy,qx)-A.omegaW*t)+A.lead/2,A.lead)-A.lead/2;
 level=max(max(abs(phase)-r30b_worm_width(max(r,A.rRoot),A),r-A.tipNew),abs(ax)-A.axialLength/2);
 core=r<=A.rRoot;coreLevel=max(r-A.rRoot,abs(ax)-A.axialLength/2);level(core)=coreLevel(core);
 along=Q.ds+abs(hc).*Q.rollStep;travel=hypot(along,Q.dz')+(hi-lo)/2;
 rbw=A.rPitch/hypot(1,A.rPitch*A.slope/A.helicalP);
 lip=1+abs(A.helicalP)./max(r-travel,A.rRoot)+abs(A.helicalP)/rbw;
 mask=level<=lip.*travel+1e-10|r-travel<=A.rRoot;
else,mask=true(Q.ns,Q.nz);end
ids=find(mask);volume=0;nCandidates=numel(ids);if isempty(ids),return;end
[i,j]=ind2sub([Q.ns,Q.nz],ids);[u,v]=ndgrid(1:Q.n,1:Q.n);u=u(:)';v=v(:)';
s=Q.s(i,u);z=Q.z(j,v);before=H(ids,:);
after=r30b_normal_update(s,z,before,t,Q.p,A,Q.opts.tooth,Q.opts.side,Q.stock);
assert(all(after(:)>=-o.normalDepthLimit_mm),'R30B:NormalLayer',...
 'Signed surface exceeded normalDepthLimit_mm; inspect geometry instead of clipping the target');
dh=before-after;after(dh<=o.heightTolerance_mm)=before(dh<=o.heightTolerance_mm);dh=before-after;
assert(all(dh(:)>=0),'R30B:Monotone','Removed material cannot be restored');
H(ids,:)=after;Q.h=reshape(H,Q.ns,Q.nz,Q.n,Q.n);
weights=Q.ds(i).*Q.dz(j)*(Q.w(u).*Q.w(v))';volume=sum(sum(weights.*dh));
if nargout>=4
 [cellIndex,combination]=find(dh>0);uq=u(:);vq=v(:);
 fi=(i(cellIndex)-1)*Q.n+uq(combination);fj=(j(cellIndex)-1)*Q.n+vq(combination);
 pick=sub2ind(size(dh),cellIndex,combination);change.index=fi+(fj-1)*(Q.ns*Q.n);
 change.before_mm=before(pick);change.after_mm=after(pick);
end
end
