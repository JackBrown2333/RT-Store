function labels=r30b_components(mask,se,ze,es,ez,lo,hi)
% Four-connected numerical supports. Diagonal contact does not merge islands.
[ns,nz]=size(mask);labels=zeros(ns,nz);seed=find(mask);queue=zeros(nnz(mask),1);reg=0;
for start=seed(:)'
 if labels(start),continue;end
 reg=reg+1;head=1;tail=1;queue(1)=start;labels(start)=reg;
 while head<=tail
  here=queue(head);head=head+1;[i,j]=ind2sub([ns,nz],here);neighbors=[i-1,j;i+1,j;i,j-1;i,j+1];
  for k=1:4
   ii=neighbors(k,1);jj=neighbors(k,2);
   if ii<1||ii>ns||jj<1||jj>nz||~mask(ii,jj)||labels(ii,jj),continue;end
   if isfinite(lo)||isfinite(hi)
    if ii~=i,s=se(max(i,ii));a=-ez*s+es*ze(j);b=-ez*s+es*ze(j+1);
    else,z=ze(max(j,jj));a=-ez*se(i)+es*z;b=-ez*se(i+1)+es*z;end
    lower=min(a,b);upper=max(a,b);
    connected=min(upper,hi)>max(lower,lo)+1e-12||...
     (abs(upper-lower)<1e-12&&lower>lo+1e-12&&lower<hi-1e-12);
    if ~connected,continue;end
   end
   labels(ii,jj)=reg;tail=tail+1;queue(tail)=sub2ind([ns,nz],ii,jj);
  end
 end
end
end
