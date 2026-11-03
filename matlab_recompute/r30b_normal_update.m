function [after,v]=r30b_normal_update(s,z,before,t,p,A,tooth,side,span)
% 固定材料点独立更新；不以h=0作为物理挡块。括定割线法失败时二分。
if nargin<7,tooth=17;end;if nargin<8,side=1;end;if nargin<9,span=.05;end
v=r30b_level(s,z,before,t,p,A,tooth,side);after=before;hit=v<0;
if ~any(hit(:)),return;end
s=s(hit);z=z(hit);b=before(hit);fh=v(hit);hi=b;lo=min(0,b-span);
fl=r30b_level(s,z,lo,t,p,A,tooth,side);width=span;
for k=1:8
 missing=fl<=0;if ~any(missing(:)),break;end
 width=2*width;lo(missing)=b(missing)-width;
 assert(all(lo(missing)>=-1),'R30:NormalLayer','法向求交超出支持的局部材料层');
 fl=r30b_level(s,z,lo,t,p,A,tooth,side);
end
assert(all(fl(:)>0),'R30:RootBracket','未能括定有符号法向交点');
active=true(size(b));answer=b;
for k=1:8
 den=fl-fh;mid=(hi+lo)/2;ok=abs(den)>1e-30;
 mid(ok)=(hi(ok).*fl(ok)-lo(ok).*fh(ok))./den(ok);
 mid=min(hi,max(lo,mid));fm=r30b_level(s,z,mid,t,p,A,tooth,side);
 done=active&abs(fm)<2e-13;answer(done)=mid(done);active(done)=false;
 inside=fm<=0;idx=active&inside;hi(idx)=mid(idx);fh(idx)=fm(idx);
 idx=active&~inside;lo(idx)=mid(idx);fl(idx)=fm(idx);
 if ~any(active(:)),break;end
end
if any(active(:))
 for k=1:34
  mid=(hi+lo)/2;fm=r30b_level(s,z,mid,t,p,A,tooth,side);inside=fm<=0;
  idx=active&inside;hi(idx)=mid(idx);idx=active&~inside;lo(idx)=mid(idx);
 end
 answer(active)=hi(active);
end
after(hit)=min(b,answer);
end
