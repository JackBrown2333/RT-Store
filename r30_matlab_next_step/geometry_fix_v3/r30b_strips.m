function rows=r30b_strips(Q,before,after,passStart,t,dt,o,capture)
if nargin<8,capture=false;end
dh=before-after;active=dh>0;rows=struct([]);if ~any(active(:)),return;end
se=Q.support_s;ze=Q.support_z;face=[Q.opts.tooth Q.opts.side];
[ii,jj]=find(active);linear=sub2ind(size(active),ii,jj);
labels=r30b_components(active,se,ze,0,0,-Inf,Inf);island=labels(linear);
area=(se(ii+1)-se(ii)).*(ze(jj+1)-ze(jj));sm=(se(ii+1)+se(ii))/2;zm=(ze(jj+1)+ze(jj))/2;
sampleS=reshape(Q.s',[],1);sampleZ=reshape(Q.z',[],1);width=o.forceStripWidth_mm;
for reg=1:max(island)
 take=find(island==reg);sc=sum(sm(take).*area(take))/sum(area(take));zc=sum(zm(take).*area(take))/sum(area(take));
 K=r30b_kinematics(Q.p,Q.A,t,face(1),face(2),sc,zc,height_at(after,se,ze,sc,zc));es=K.es;ez=K.ez;
 uu=[-ez*se(ii(take))+es*ze(jj(take)),-ez*se(ii(take)+1)+es*ze(jj(take)),...
  -ez*se(ii(take)+1)+es*ze(jj(take)+1),-ez*se(ii(take))+es*ze(jj(take)+1)];
 umin=min(uu,[],2);umax=max(uu,[],2);lower=floor(min(umin)/width)*width;upper=ceil(max(umax)/width)*width;
 for ulo=lower:width:upper-width/2
  candidates=take(umax>ulo&umin<ulo+width);if isempty(candidates),continue;end
  slab=false(size(active));slab(linear(candidates))=true;
  lab=r30b_components(slab,se,ze,es,ez,ulo,ulo+width);subids=lab(linear(candidates));
  for part=1:max(subids)
   ids=candidates(subids==part);i=ii(ids);j=jj(ids);lin=linear(ids);count=numel(ids);V=zeros(count,4,11);
   V(:,:,1)=[se(i),se(i+1),se(i+1),se(i)];V(:,:,2)=[ze(j),ze(j),ze(j+1),ze(j+1)];
   for col=[3 6],V(:,:,col)=repmat(dh(lin),1,4);end
   V(:,:,4)=repmat(before(lin),1,4);V(:,:,5)=repmat(after(lin),1,4);
   V(:,:,7)=repmat(passStart(lin),1,4);V(:,:,8)=repmat(passStart(lin)-after(lin),1,4);
   V(:,:,9)=repmat(sampleS(i),1,4);V(:,:,10)=repmat(sampleZ(j),1,4);V(:,:,11)=-ez*V(:,:,1)+es*V(:,:,2);
   [B,bc]=r30b_clip(V,ones(count,1)*4,11,ulo,1);[B,bc]=r30b_clip(B,bc,11,ulo+width,-1);
   [ar,cs,cz,vol]=r30b_geometry(B);keep=ar>0;B=B(keep,:,:);bc=bc(keep);ar=ar(keep);cs=cs(keep);cz=cz(keep);vol=vol(keep);
   if isempty(ar),continue;end
   totalArea=sum(ar);ss=sum(ar.*cs)/totalArea;zz=sum(ar.*cz)/totalArea;hh=height_at(after,se,ze,ss,zz);
   center=r30b_kinematics(Q.p,Q.A,t,face(1),face(2),ss,zz,hh);
   low=min(min(B(:,:,11)));high=max(max(B(:,:,11)));b=high-low;assert(b>0);
   [positions,normals,tangents]=r30b_kinematics_many(cs,cz,mean(B(:,:,5),2),t,Q.p,Q.A,face(1),face(2));
   weights=ar/totalArea;unitN=-sum(normals.*weights,1);unitT=sum(tangents.*weights,1);
   momentN=sum(cross(positions,-normals,2).*weights,1);momentT=sum(cross(positions,tangents,2).*weights,1);
   rec=struct('tooth_id',face(1),'flank_side',face(2),'region_id',reg,'strip_component',part,...
    't_s',t,'dt_s',dt,'b_mm',b,'lg_mm',totalArea/b,'contact_area_mm2',totalArea,'dV_mm3',sum(vol),...
    'ae_pass_mm',max(max(B(:,:,8))),'dh_step_mm',max(max(B(:,:,6))),...
    'stock_nominal_mm',Q.stock,'vc_mmps',center.vc,'contact_s_mm',ss,'contact_z_mm',zz,...
    'x_mm',center.point(1),'y_mm',center.point(2),'z_mm',center.point(3),...
    'normal_response',unitN,'tangent_response',unitT,'normal_moment_response',momentN,'tangent_moment_response',momentT);
   if capture
    rec.geometry=annotations(B,bc,center,Q,t,ss,zz,es,ez,low,high,rec.lg_mm,o);
   end
   if isempty(rows),rows=rec;else,rows(end+1,1)=rec;end
  end
 end
end
end
function h=height_at(H,se,ze,s,z)
if ~isscalar(s) || ~isscalar(z) || ~isfinite(s) || ~isfinite(z)
 error('R30B:StripCenter','Strip center coordinates must be finite scalars.');
end
scale=max([1;abs(se(:));abs(ze(:))]);
tol=128*eps(scale);
if s<se(1)-tol || s>se(end)+tol || z<ze(1)-tol || z>ze(end)+tol
 error('R30B:StripCenterDomain','Strip center is outside the material support.');
end
s=min(max(s,se(1)),se(end));
z=min(max(z,ze(1)),ze(end));
i=find(se<=s,1,'last');
j=find(ze<=z,1,'last');
if isempty(i) || isempty(j)
 error('R30B:StripCenterIndex','Material support lookup returned an empty index.');
end
i=min(i,size(H,1));
j=min(j,size(H,2));
h=H(i,j);
if ~isscalar(h) || ~isfinite(h)
 error('R30B:StripCenterHeight','Strip center height must be a finite scalar.');
end
end
function g=annotations(B,bc,K,Q,t,s,z,es,ez,lo,hi,lg,o)
[~,~,~,~,tx,ty]=r30b_surface(s,Q.opts.tooth,Q.opts.side,Q.p,Q.A);th=Q.A.phi0+Q.A.omegaG*t;
ts=[tx*cos(th)-ty*sin(th);tx*sin(th)+ty*cos(th);0];along=es*ts+ez*[0;0;1];across=-ez*ts+es*[0;0;1];
uc=-ez*s+es*z;point=K.point(:)';
g=struct('lengthCoordinates','target-surface developed s,z; thin-layer measure',...
 'lgDefinition','area_over_transverse_width_equivalent_mean_chord',...
 'b_segment_world_mm',[point+(lo-uc)*across';point+(hi-uc)*across'],...
 'lg_equivalent_segment_world_mm',[point-lg/2*along';point+lg/2*along'],...
 'stripDirection_sz',[es ez],'transverseRange_mm',[lo hi]);
segments=r30b_section(B,es,ez,(lo+hi)/2);g.actualChordSegments_v_mm=segments;
g.actualChordLine_u_mm=(lo+hi)/2;
g.actualChordSegments_sz=zeros(size(segments,1),4);
for k=1:size(segments,1)
 v=segments(k,:);u=(lo+hi)/2;g.actualChordSegments_sz(k,:)=[es*v(1)-ez*u,ez*v(1)+es*u,es*v(2)-ez*u,ez*v(2)+es*u];
end
[~,idx]=max(B(:,1,8));sp=B(idx,1,9);zp=B(idx,1,10);
a=r30b_kinematics(Q.p,Q.A,t,Q.opts.tooth,Q.opts.side,sp,zp,B(idx,1,7));
b=r30b_kinematics(Q.p,Q.A,t,Q.opts.tooth,Q.opts.side,sp,zp,B(idx,1,5));
g.ae_pass_anchor_sz=[sp zp];g.ae_pass_segment_world_mm=[a.point';b.point'];
[~,idx]=max(B(:,1,6));sp=B(idx,1,9);zp=B(idx,1,10);
a=r30b_kinematics(Q.p,Q.A,t,Q.opts.tooth,Q.opts.side,sp,zp,B(idx,1,4));
b=r30b_kinematics(Q.p,Q.A,t,Q.opts.tooth,Q.opts.side,sp,zp,B(idx,1,5));
g.dh_step_segment_world_mm=[a.point';b.point'];
if o.saveStripPolygons,g.clippedSupportPolygons=B(:,:,1:8);g.vertexCounts=bc;end
end
