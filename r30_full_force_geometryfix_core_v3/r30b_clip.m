function [out,cnt]=r30b_clip(Q,counts,col,bound,sense)
% Batched convex-polygon clipping. Vertices are never capped at seven.
[n,m,d]=size(Q);if n==0,out=Q;cnt=counts;return;end
cap=m+1;out=zeros(n,cap,d);cnt=zeros(n,1);rows=(1:n)';
for k=1:m
 prev=k-1+zeros(n,1);if k==1,prev=counts;end
 P=zeros(n,d);N=zeros(n,d);
 for f=1:d
  a=Q(:,:,f);P(:,f)=a(sub2ind([n,m],rows,prev));N(:,f)=a(:,k);
 end
 vp=sense*(P(:,col)-bound);vn=sense*(N(:,col)-bound);valid=k<=counts;
 cross=valid&xor(vp>0,vn>0);den=vp-vn;frac=zeros(n,1);ok=abs(den)>1e-30;frac(ok)=vp(ok)./den(ok);
 I=P+(N-P).*frac;rr=find(cross);ix=sub2ind([n,cap],rr,cnt(rr)+1);
 for f=1:d,a=out(:,:,f);a(ix)=I(rr,f);out(:,:,f)=a;end
 cnt(rr)=cnt(rr)+1;
 rr=find(valid&vn>0);ix=sub2ind([n,cap],rr,cnt(rr)+1);
 for f=1:d,a=out(:,:,f);a(ix)=N(rr,f);out(:,:,f)=a;end
 cnt(rr)=cnt(rr)+1;
end
good=cnt>=3;out=out(good,:,:);cnt=cnt(good);n=sum(good);
if n
 rr=repmat((1:n)',1,cap);cc=min(repmat(1:cap,n,1),cnt);
 for f=1:d,a=out(:,:,f);out(:,:,f)=reshape(a(sub2ind([n,cap],rr,cc)),n,cap);end
end
end
