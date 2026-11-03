function segments=r30b_section(B,es,ez,u0)
% Actual intersection of a line u=u0 with the numerical-support polygon union.
U=-ez*B(:,:,1)+es*B(:,:,2);V=es*B(:,:,1)+ez*B(:,:,2);
Un=U(:,[2:end 1]);Vn=V(:,[2:end 1]);den=Un-U;
crossing=(min(U,Un)<=u0+1e-12)&(max(U,Un)>=u0-1e-12)&abs(den)>1e-14;
vv=V+(u0-U)./den.*(Vn-V);vv(~crossing)=NaN;
on=abs(U-u0)<1e-12;vertices=V;vertices(~on)=NaN;vv=[vv,vertices];
low=vv;low(isnan(low))=Inf;high=vv;high(isnan(high))=-Inf;
intervals=[min(low,[],2),max(high,[],2)];intervals=intervals(intervals(:,2)>intervals(:,1)+1e-12,:);
intervals=sortrows(intervals,1);segments=zeros(0,2);
for k=1:size(intervals,1)
 if isempty(segments)||intervals(k,1)>segments(end,2)+1e-10,segments(end+1,:)=intervals(k,:); %#ok<AGROW>
 else,segments(end,2)=max(segments(end,2),intervals(k,2));end
end
end
