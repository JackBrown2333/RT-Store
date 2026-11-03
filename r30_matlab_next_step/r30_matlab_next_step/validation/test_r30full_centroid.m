function report=test_r30full_centroid(originalCore,fixedCore,sourceM)
% Slivers, translation, clipped supports, and full force-chain height lookup.
addpath(fixedCore,'-begin');clear r30b_geometry r30b_strips;
assert(strcmp(fileparts(which('r30b_geometry')),fixedCore));
Q=zeros(1,4,6);Q(1,:,1)=[0 1e-4 1e-4 0];Q(1,:,2)=[-10 -10 -10+1e-12 -10+1e-12];Q(:,:,6)=.01;
[ar,s,z,v]=r30b_geometry(Q);[oldArea,oldS,oldZ]=old_geometry(Q);
assert(oldArea>0&&oldZ<-10&&s>=0&&s<=1e-4&&z>=-10&&z<=-10+1e-12);
assert(abs(s-5e-5)<1e-18&&abs(z-(-10+(.5e-12)))<4e-15&&abs(v-.01*ar)<1e-30);
% A translated ordinary polygon must have identical area and translated center.
T=zeros(1,4,6);T(1,:,1)=[0 .125 .125 0];T(1,:,2)=[0 0 .25 .25];T(:,:,6)=.02;
[a0,s0,z0,v0]=r30b_geometry(T);T(:,:,1)=T(:,:,1)+128;T(:,:,2)=T(:,:,2)-64;
[a1,s1,z1,v1]=r30b_geometry(T);assert(a0==a1&&v0==v1&&abs(s1-s0-128)<1e-14&&abs(z1-z0+64)<1e-14);
% Single-support and face-edge force paths, including a signed negative H.
D=load(sourceM,'M');P=D.M.parameters;clear D;
H=ones(2,2)*.05;se=[0;1e-4;2e-4];ze=[-10;-9.999;-9.998];
i=find(se<=oldS,1,'last');j=find(ze<=oldZ,1,'last');emptyHeight=H(i,j);assert(isempty(emptyHeight));
oldEmptyHeightError=false;oldFailureIdentifier='';oldFailureMessage='';
try,r30b_kinematics(P.p,P.installation,2.25734,1,1,oldS,oldZ,emptyHeight);
catch e,oldEmptyHeightError=true;oldFailureIdentifier=e.identifier;oldFailureMessage=e.message;end
assert(oldEmptyHeightError,'R30FULL:Regression','未复现空高度导致的运动学计算错误。');
i=find(se<=s,1,'last');j=find(ze<=z,1,'last');K=r30b_kinematics(P.p,P.installation,2.25734,1,1,s,z,H(i,j));
assert(numel(K.point)==3&&all(isfinite(K.point)));
co=r30b_options(struct('ns',12,'nz',24,'order',2,'forceStripWidth_mm',.05,'plot',false));
q=r30b_init(P,co,[1 1]);after=r30b_heights(q);before=after;before(1,1)=before(1,1)+1e-6;
rows=r30b_force(r30b_strips(q,before,after,before,2.25734,co.dt_s,co,false),co.coefficients);
assert(~isempty(rows)&&all(isfinite([rows.Fn_N])));
for signZ=[1 size(after,2)]
 after(:)=.05;after(1,signZ)=-.001;before=after;before(1,signZ)=before(1,signZ)+1e-6;
 rows=r30b_force(r30b_strips(q,before,after,before,2.25734,co.dt_s,co,false),co.coefficients);
 assert(~isempty(rows)&&all(isfinite([rows.Fn_N]))&&after(1,signZ)==-.001);
end
report=struct('oldSliverCentroid',[oldS oldZ],'fixedSliverCentroid',[s z],...
 'sliverArea_mm2',ar,'reproducesOldOutOfBoundsLookup',true,'localCoordinateCentroidInside',true,...
 'translationInvariantOrdinaryPolygon',true,'edgeForceChainsFinite',true,'negativeMaterialHeightPreserved',true,...
 'originalEmptyHeightFailureReproduced',oldEmptyHeightError,'originalFailureIdentifier',oldFailureIdentifier,...
 'originalFailureMessage',oldFailureMessage,'fixedKinematicsHasThreeCoordinates',true,...
 'actualUserFailureSnapshotAvailable',false,'originalCoreDir',originalCore);
disp(report);
end
function [area,s,z]=old_geometry(Q)
x=Q(:,:,1);y=Q(:,:,2);xn=x(:,[2:end 1]);yn=y(:,[2:end 1]);c=x.*yn-xn.*y;d=sum(c,2);area=abs(d)/2;
s=sum((x+xn).*c,2)./(3*d);z=sum((y+yn).*c,2)./(3*d);
end
