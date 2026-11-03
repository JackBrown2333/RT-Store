function Q=r30b_quadrature(P,opts)
% 固定材料Gauss积分点。权重支撑矩形互不重叠，其面积等于积分权重。
if nargin<2,opts=struct();end
defaults=struct('order',2,'ns',400,'nz',800,'z0',-10,'z1',10,'tooth',17,'side',1);
names=fieldnames(defaults);
for j=1:numel(names),if ~isfield(opts,names{j}),opts.(names{j})=defaults.(names{j});end;end
n=opts.order;
if n==1,q=.5;w=1;
else
 k=(1:n-1)';beta=k./sqrt(4*k.^2-1);[V,D]=eig(diag(beta,1)+diag(beta,-1));
 [q,ix]=sort(diag(D));w=(V(1,ix).^2)';q=(q+1)/2;
end
p=P.p;A=P.installation;A.axisW=A.axisW(:);A.axisY=A.axisY(:);
rb=p.calc.gear.baseDiameter_mm/2;ra=p.calc.gear.tipDiameter_mm/2;smax=(ra^2-rb^2)/(2*rb);
se=smax*linspace(0,1,opts.ns+1).^2;ze=linspace(opts.z0,opts.z1,opts.nz+1);
ds=diff(se)';dz=diff(ze)';s=se(1:end-1)'+ds*q';z=ze(1:end-1)'+dz*q';
cum=[0;cumsum(w)];ss=se(1:end-1)'+ds*cum(1:end-1)';zz=ze(1:end-1)'+dz*cum(1:end-1)';
support_s=[reshape(ss',[],1);se(end)];support_z=[reshape(zz',[],1);ze(end)];
Q=struct('p',p,'A',A,'opts',opts,'n',n,'ns',opts.ns,'nz',opts.nz,...
 'rb',rb,'smax',smax,'se',se,'ze',ze,'ds',ds,'dz',dz,'q',q,'w',w,...
 's',s,'z',z,'support_s',support_s,'support_z',support_z,...
 'h',ones(opts.ns,opts.nz,n,n)*P.cfg.stock_um/1000,'stock',P.cfg.stock_um/1000);
end
