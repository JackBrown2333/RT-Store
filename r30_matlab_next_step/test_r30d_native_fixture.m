function report=test_r30d_native_fixture(sourceM,coreDir,fixtureDir)
% Synthetic short-window fixture for the diagnostic reader; NOT project data.
addpath(coreDir);D=load(sourceM,'M');parameters=D.M.parameters;clear D;
o=r30b_options(struct('ns',24,'nz',48,'order',2,'dt_s',1.25e-5,'substeps',2,'plot',false));
nt=3;t0=15.0139;faces=[17 1;17 2];initialState=struct('t_s',t0,'faces',{{}});Qs=cell(2,1);
for side=1:2
 Q=r30b_init(parameters,o,faces(side,:));Q.passStart=Q.h;Q.inPass=true;
 [Q,~]=r30b_step(Q,t0,o);Qs{side}=Q;H=r30b_heights(Q);temp=Q;temp.h=Q.passStart;pass=r30b_heights(temp);
 initialState.faces{side}=struct('tooth_id',17,'flank_side',side,'height_mm',H,...
  'passStartHeight_mm',pass,'support_s_mm',Q.support_s,'support_z_mm',Q.support_z,...
  'sample_s_mm',reshape(Q.s',[],1),'sample_z_mm',reshape(Q.z',[],1),'stock_nominal_mm',Q.stock);
end
leaf=fullfile(fixtureDir,'animation_export','replay','tooth_17');mkdir(fullfile(leaf,'frames'));
save(fullfile(leaf,'initial_state.mat'),'initialState','-v7');
M=struct('parameters',parameters,'toothIds',17,'faces',faces,'calculatedHistoryForTooth',true,...
 'forceDt_s',o.dt_s,'time_s',t0+(1:nt)'*o.dt_s,'window_s',[t0 t0+nt*o.dt_s],...
 'force_N',zeros(nt,3),'moment_Nmm',zeros(nt,3),'F_resultant_N',zeros(nt,1),...
 'toothForce_N',zeros(nt,1,3),'toothMoment_Nmm',zeros(nt,1,3),...
 'perToothResultFiles',{{fullfile('tooth_17','r30b_result.mat')}},'coefficients',o.coefficients);
R30B=struct('faces',faces,'parameters',parameters,'options',o,'forceCoefficients',o.coefficients,...
 'momentOrigin_mm',[0 0 0],'initialStateFile','initial_state.mat','frameIndices',1:nt,'frameFiles',{{}},'facePoses',{{}});
for k=1:nt
 ta=t0+(k-1)*o.dt_s;t=t0+k*o.dt_s;frame=struct('index',k,'t_start_s',ta,'t_s',t,'dt_s',o.dt_s,'faces',{{}});totalF=[0 0 0];totalM=totalF;
 for side=1:2
  Q=Qs{side};before=r30b_heights(Q);changes=cell(o.substeps,1);
  for sub=1:o.substeps,[Q,~,~,changes{sub}]=r30b_step(Q,ta+sub*o.dt_s/o.substeps,o);end
  after=r30b_heights(Q);temp=Q;temp.h=Q.passStart;pass=r30b_heights(temp);dh=before-after;ids=find(dh>0);
  rows=r30b_force(r30b_strips(Q,before,after,pass,t,o.dt_s,o,true),o.coefficients);
  area=sum(Q.weightsFine(dh>0));volume=sum(Q.weightsFine(:).*dh(:));sumrow=r30b_sum(rows,o.referenceOutward,t,o.dt_s,volume,area);
  R30B.facePoses{side}(k)=sumrow;
  frame.faces{side}=struct('stateUpdates',struct('since_t_s',ta,'index',ids,'height_mm',after(ids)),...
   'removal',struct('index',ids,'before_mm',before(ids),'after_mm',after(ids),'passStart_mm',pass(ids)),...
   'strips',rows,'substepUpdates',{changes});
  totalF=totalF+sumrow.force_N;totalM=totalM+sumrow.moment_Nmm;Qs{side}=Q;
 end
 R30B.poses(k)=struct('t_s',t,'force_N',totalF,'moment_Nmm',totalM);
 M.force_N(k,:)=totalF;M.moment_Nmm(k,:)=totalM;M.F_resultant_N(k)=norm(totalF);
 M.toothForce_N(k,1,:)=reshape(totalF,1,1,3);M.toothMoment_Nmm(k,1,:)=reshape(totalM,1,1,3);
 name=sprintf('frames/frame_%06d.mat',k);save(fullfile(leaf,name),'frame','-v7');R30B.frameFiles{k}=name;
end
save(fullfile(leaf,'r30b_result.mat'),'R30B','-v7');save(fullfile(fixtureDir,'multitooth_pilot.mat'),'M','-v7');
files=dir(fullfile(coreDir,'*.m'));signature='';for j=1:numel(files),signature=[signature files(j).name char(10) fileread(fullfile(coreDir,files(j).name))];end
S=struct('phase','done','meta',struct('codeSignature',signature));mkdir(fullfile(fixtureDir,'tooth_17'));save(fullfile(fixtureDir,'tooth_17','r30b_checkpoint.mat'),'S','-v7');
A=run_r30_contact_direction_audit(fixtureDir,struct('coreDir',coreDir,'sourceIndices',1:nt,'outputDir',fullfile(fixtureDir,'audit_result')));
% Independent finite differences at actual force-quadrature coordinates.
data=dlmread(fullfile(A.outputDir,'force_quadrature_points.csv'),',',1,0);p=parameters.p;a=parameters.installation;
take=1:min(40,size(data,1));d=data(take,:);P=zeros(numel(take),3);
for j=1:numel(take),P(j,:)=r30b_kinematics_many(d(j,7),d(j,8),d(j,9),d(j,2),p,a,d(j,3),d(j,4));end
G=r30d_wheel_geometry(P,d(:,2),a);epsX=1e-5;fd=zeros(size(P));
for c=1:3,plus=P;minus=P;plus(:,c)=plus(:,c)+epsX;minus(:,c)=minus(:,c)-epsX;gp=r30d_wheel_geometry(plus,d(:,2),a);gm=r30d_wheel_geometry(minus,d(:,2),a);fd(:,c)=(gp.level-gm.level)/(2*epsX);end
gap=G.branchGap>1e-3;assert(any(gap));err=max(abs(reshape(fd(gap,:)-G.gradient(gap,:),[],1)));assert(err<1e-5);
epsT=1e-7;gp=r30d_wheel_geometry(P,d(:,2)+epsT,a);gm=r30d_wheel_geometry(P,d(:,2)-epsT,a);errT=max(abs((gp.level(gap)-gm.level(gap))/(2*epsT)-G.partialTimeDerivative(gap)));assert(errT<1e-3);
report=struct('synthetic_fixture_only',true,'full_native_reader_and_reconstruction',true,...
 'forceQuadraturePoints',A.forceQuadraturePointCount,'materialRootEvents',A.materialRootEventCount,...
 'spatial_gradient_finite_difference_error',err,'time_derivative_finite_difference_error',errT,...
 'max_errors',A.maxErrors,'actual_user_full_resolution_diagnostic_run',false);
fid=fopen(fullfile(fixtureDir,'test_report.json'),'w');fprintf(fid,'%s',jsonencode(report));fclose(fid);disp(report);
end
