function Aout=run_r30_contact_direction_audit(resultDir,user)
% RUN_R30_CONTACT_DIRECTION_AUDIT Read native animation replay frames.
% Put this package beside matlab_recompute/ and the full multitooth result.
%   Aout=run_r30_contact_direction_audit;
% Reads existing animation_export/replay and checks force-quadrature versus
% actual fixed-material root points at their respective substep times.
% Does not change source results, core code, force definitions or coefficients.
root=fileparts(mfilename('fullpath'));projectRoot=root;
if ~isfolder(fullfile(root,'matlab_recompute'))&&isfolder(fullfile(fileparts(root),'matlab_recompute')),projectRoot=fileparts(root);end
if nargin<1||isempty(resultDir),resultDir=fullfile(projectRoot,'r30_multitooth_pilot_15_v2_160x320');end
if nargin<2,user=struct();end
o=struct('coreDir',fullfile(projectRoot,'matlab_recompute'),...
 'replayDir','','outputDir','','sourceIndices',[],'progressEvery_s',5);
names=fieldnames(user);for j=1:numel(names),assert(isfield(o,names{j}),'R30D:Option','Unknown option: %s',names{j});o.(names{j})=user.(names{j});end
assert(isfolder(o.coreDir),'R30D:Core','请把脚本放在matlab_recompute同级，或用coreDir指定原核心目录。');
addpath(o.coreDir);D=load(fullfile(resultDir,'multitooth_pilot.mat'),'M');M=D.M;clear D;
if isempty(o.replayDir),o.replayDir=fullfile(resultDir,'animation_export','replay');end
assert(isfolder(o.replayDir),'R30D:Replay','缺少原生重放目录：%s。需要本机animation_export/replay，不是精简动画ZIP。',o.replayDir);
signature=core_signature(o.coreDir);ng=numel(M.toothIds);leaves=cell(ng,1);available=[];
for g=1:ng
 if ~M.calculatedHistoryForTooth(g),continue;end
 leaf=fullfile(o.replayDir,sprintf('tooth_%02d',M.toothIds(g)));rp=fullfile(leaf,'r30b_result.mat');
 assert(isfile(rp),'R30D:Replay','缺少已完成重放结果：%s',rp);D=load(rp,'R30B');R=D.R30B;clear D;
 validate_leaf(R,M,g);leaves{g}=struct('directory',leaf,'R',R);
 if isempty(available),available=double(R.frameIndices(:)');else,available=intersect(available,double(R.frameIndices(:)'));end
 sourceRelative=strrep(M.perToothResultFiles{g},'\',filesep);cp=fullfile(fileparts(fullfile(resultDir,sourceRelative)),'r30b_checkpoint.mat');
 assert(isfile(cp),'R30D:Checkpoint','缺少原逐齿检查点，无法核对核心代码：%s',cp);D=load(cp,'S');
 assert(strcmp(D.S.phase,'done')&&strcmp(D.S.meta.codeSignature,signature),'R30D:Code','核心代码与原计算不同，请恢复原matlab_recompute后诊断。');clear D;
end
assert(~isempty(available),'R30D:Frames','没有各齿共同保存的原生帧。');
tc=M.moment_Nmm(:,3);[~,kp]=max(tc);[~,kn]=min(tc);[~,kf]=max(M.F_resultant_N);
crossings=find(tc(2:end).*tc(1:end-1)<0)+1;
if isempty(o.sourceIndices),wanted=unique([kp kn kf crossings(1:min(3,numel(crossings)))']);else,wanted=unique(o.sourceIndices(:)');end
assert(all(wanted>=1&wanted<=numel(M.time_s)&mod(wanted,1)==0));
actual=zeros(size(wanted));for j=1:numel(wanted),[~,q]=min(abs(available-wanted(j)));actual(j)=available(q);end
selected=unique(actual);selection=[wanted(:),actual(:),1e6*(M.time_s(actual(:))-M.time_s(wanted(:)))];
if isempty(o.outputDir),o.outputDir=unique_folder(fullfile(resultDir,['contact_direction_audit_' datestr(now,'yyyymmdd_HHMMSS')]));end
assert(~isfolder(o.outputDir),'R30D:Output','输出目录已存在，请换新目录。');mkdir(o.outputDir);[ok,attr]=fileattrib(o.outputDir);assert(ok);out=attr.Name;
fprintf('方向诊断：%d个代表帧；只读取现有原生重放帧。\n',numel(selected));
forcePoints=zeros(0,22);rootPoints=zeros(0,16);frameTotals=zeros(numel(selected),9);
frameTotals(:,1)=selected;frameTotals(:,2)=M.time_s(selected);
maxerr=struct('force_N',0,'moment_Nmm',0,'normal_response',0,'tangent_response',0,...
 'normal_moment_response_mm',0,'tangent_moment_response_mm',0,'root_level',0,'rigid_transport',0,'force_point_level_match',0,'frame_volume_mm3',0);
clock=tic;
for g=1:ng
 if isempty(leaves{g}),continue;end
 R=leaves{g}.R;leaf=leaves{g}.directory;D=load(fullfile(leaf,R.initialStateFile),'initialState');initial=D.initialState;clear D;
 nf=size(R.faces,1);H=cell(nf,1);Qall=H;
 for side=1:nf
  q=initial.faces{side};Q=r30b_init(R.parameters,R.options,R.faces(side,:));
  assert(isequal(size(q.height_mm),[Q.ns*Q.n Q.nz*Q.n])&&q.tooth_id==R.faces(side,1)&&q.flank_side==R.faces(side,2));
  assert(max(abs(q.support_s_mm(:)-Q.support_s(:)))<1e-13&&max(abs(q.support_z_mm(:)-Q.support_z(:)))<1e-13);
  H{side}=q.height_mm;Qall{side}=Q;
 end
 lastTime=initial.t_s;
 for j=1:numel(R.frameFiles)
  k=double(R.frameIndices(j));if k>max(selected),break;end
  D=load(fullfile(leaf,strrep(R.frameFiles{j},'\',filesep)),'frame');fr=D.frame;clear D;
  assert(fr.index==k&&abs(fr.t_s-M.time_s(k))<1e-10);isSelected=ismember(k,selected);
  for side=1:nf
   x=fr.faces{side};u=x.stateUpdates;assert(abs(u.since_t_s-lastTime)<1e-10);
   ids=double(u.index(:));values=u.height_mm(:);old=H{side};
   assert(numel(ids)==numel(values)&&numel(unique(ids))==numel(ids)&&all(ids>=1&ids<=numel(old))&&all(isfinite(values)));
   assert(all(values<=old(ids)+1e-13),'R30D:Replay','材料增量恢复了已去除材料。');old(ids)=values;H{side}=old;
   if ~isSelected,continue;end
   Q=Qall{side};after=old;before=after;pass=after;remove=x.removal;rid=double(remove.index(:));
   assert(max0(abs(after(rid)-remove.after_mm(:)))<1e-12);before(rid)=remove.before_mm(:);pass(rid)=remove.passStart_mm(:);
   assert(all(before(:)>=after(:)));Q.h=unpack(after,Q);capture=R.options;capture.saveStripPolygons=true;
   rows=r30b_force(r30b_strips(Q,before,after,pass,fr.t_s,fr.dt_s,capture,true),capture.coefficients);
   stored=x.strips;assert(numel(rows)==numel(stored),'R30D:Strips','重建条带数量不同。');
   sumF=[0 0 0];sumM=sumF;pn=0;pt=0;
   for r=1:numel(rows)
    rr=rows(r);ss=stored(r);F=[rr.Fx_N rr.Fy_N rr.Fz_N];mm=[rr.Mx_Nmm rr.My_Nmm rr.Mz_Nmm];
    maxerr.force_N=max(maxerr.force_N,max(abs(F-[ss.Fx_N ss.Fy_N ss.Fz_N])));
    maxerr.moment_Nmm=max(maxerr.moment_Nmm,max(abs(mm-[ss.Mx_Nmm ss.My_Nmm ss.Mz_Nmm])));
    assert(maxerr.force_N<1e-7&&maxerr.moment_Nmm<1e-5,'R30D:Reconstruct','原力或力矩不能重现，已停止诊断。');
    B=rr.geometry.clippedSupportPolygons;[area,s,z]=r30b_geometry(B);h=mean(B(:,:,5),2);weight=area/sum(area);
    [point,N,E]=r30b_kinematics_many(s,z,h,fr.t_s,Q.p,Q.A,rr.tooth_id,rr.flank_side);G=r30d_wheel_geometry(point,fr.t_s,Q.A);
    rel=relative_velocity(point,Q.A);inward=-N;dFn=rr.Fn_N*weight;dFt=rr.Ft_N*weight;
    nrel=sum(inward.*rel,2);wrel=sum(G.outwardNormal.*rel,2);normalPower=-dFn.*nrel/1000;tangentPower=-dFt.*sum(E.*rel,2)/1000;
    angle=acosd(max(-1,min(1,sum(inward.*G.outwardNormal,2))));
    nresp=sum(bsxfun(@times,inward,weight),1);tresp=sum(bsxfun(@times,E,weight),1);
    nmresp=sum(bsxfun(@times,cross(point,inward,2),weight),1);tmresp=sum(bsxfun(@times,cross(point,E,2),weight),1);
    for pair={{'normal_response',nresp},{'tangent_response',tresp},{'normal_moment_response_mm',nmresp},{'tangent_moment_response_mm',tmresp}}
     v=pair{1};fieldName=strrep(v{1},'_mm','');maxerr.(v{1})=max(maxerr.(v{1}),max(abs(v{2}-rr.(fieldName))));
    end
    level=r30b_level(s,z,h,fr.t_s,Q.p,Q.A,rr.tooth_id,rr.flank_side);
    maxerr.force_point_level_match=max(maxerr.force_point_level_match,max0(abs(level-G.level)));
    maxerr.rigid_transport=max(maxerr.rigid_transport,max0(abs(G.rigidTransportResidual)));
    count=numel(area);c=@(v)repmat(v,count,1);
    forcePoints=[forcePoints;[c(k),c(fr.t_s),c(rr.tooth_id),c(rr.flank_side),c(r),(1:count)',s,z,h,area,dFn,dFt,angle,G.level,G.branch,G.branchGap,nrel,wrel,normalPower,tangentPower,normalPower+tangentPower,point(:,3)]]; %#ok<AGROW>
    sumF=sumF+F;sumM=sumM+mm;pn=pn+sum(normalPower);pt=pt+sum(tangentPower);
   end
   target=R.facePoses{side}(k);maxerr.force_N=max(maxerr.force_N,max(abs(sumF-target.force_N)));maxerr.moment_Nmm=max(maxerr.moment_Nmm,max(abs(sumM-target.moment_Nmm)));
   volume=sum(Q.weightsFine(:).*(before(:)-after(:)));maxerr.frame_volume_mm3=max(maxerr.frame_volume_mm3,abs(volume-target.dV_mm3));
   pos=find(selected==k);frameTotals(pos,3:9)=frameTotals(pos,3:9)+[pn pt pn+pt sumF sumM(3)/1000];
   % Original material-column cuts: correct coordinates and each cut's actual substep time.
   assert(isfield(x,'substepUpdates'),'R30D:Substeps','帧缺少原子步更新，不能冒充真实求交点诊断。');
   for sub=1:numel(x.substepUpdates)
    change=x.substepUpdates{sub};if isempty(change),continue;end
    si=double(change.index(:));if isempty(si),continue;end
    [ii,jj]=ind2sub(size(after),si);sampleS=reshape(Q.s',[],1);sampleZ=reshape(Q.z',[],1);
    s=sampleS(ii);z=sampleZ(jj);h=change.after_mm(:);ts=change.t_s;
    [point,N]=r30b_kinematics_many(s,z,h,ts,Q.p,Q.A,R.faces(side,1),R.faces(side,2));G=r30d_wheel_geometry(point,ts,Q.A);rel=relative_velocity(point,Q.A);
    coreLevel=r30b_level(s,z,h,ts,Q.p,Q.A,R.faces(side,1),R.faces(side,2));maxerr.root_level=max(maxerr.root_level,max0(abs(coreLevel)));
    maxerr.rigid_transport=max(maxerr.rigid_transport,max0(abs(G.rigidTransportResidual)));
    angle=acosd(max(-1,min(1,sum((-N).*G.outwardNormal,2))));nn=numel(si);c=@(v)repmat(v,nn,1);
    rootPoints=[rootPoints;[c(k),c(fr.t_s),c(ts),c(R.faces(side,1)),c(R.faces(side,2)),si,s,z,h,change.before_mm(:)-h,G.level,G.branch,G.branchGap,angle,sum((-N).*rel,2),sum(G.outwardNormal.*rel,2)]]; %#ok<AGROW>
   end
  end
  lastTime=fr.t_s;
  if toc(clock)>o.progressEvery_s,fprintf('  第%d齿：已读取帧%d；积分点%d，子步求交点%d。\n',M.toothIds(g),k,size(forcePoints,1),size(rootPoints,1));clock=tic;end
 end
end
assert(~isempty(forcePoints)&&~isempty(rootPoints),'R30D:Empty','代表帧没有有效切削点。');
expectedPower=zeros(numel(selected),1);
for j=1:numel(selected)
 k=selected(j);F=M.force_N(k,:);mm=M.moment_Nmm(k,:);a=M.parameters.installation;
 expectedPower(j)=(a.omegaG*mm(3)-a.feed*F(3)+dot(a.omegaW*a.axisW,-mm+cross([a.center 0 0],F)))/1000;
end
maxerr.total_power_W=max(abs(frameTotals(:,5)-expectedPower));maxerr.global_force_N=max0(abs(frameTotals(:,6:8)-M.force_N(selected,:)));
maxerr.global_C_torque_Nm=max(abs(frameTotals(:,9)-M.moment_Nmm(selected,3)/1000));
assert(maxerr.force_N<1e-7&&maxerr.moment_Nmm<1e-5&&maxerr.frame_volume_mm3<1e-9&&maxerr.total_power_W<1e-6&&maxerr.global_force_N<1e-7);
assert(maxerr.force_point_level_match<1e-8&&maxerr.rigid_transport<1e-7&&maxerr.root_level<1e-7,'R30D:Geometry','隐式几何或真实求交点检查不符，请保留输出目录供排查。');
write_csv(fullfile(out,'selected_frames.csv'),{'wanted_source_index','used_source_index','time_offset_us'},selection);
write_csv(fullfile(out,'force_quadrature_points.csv'),{'source_index','time_s','tooth','side','strip','polygon','s_mm','z_mm','h_mm','area_mm2','Fn_share_N','Ft_share_N','normal_angle_deg','wheel_implicit_level_NOT_distance','branch','branch_gap','target_inward_dot_vrel_mmps','wheel_outward_dot_vrel_mmps','normal_pair_load_power_W','tangent_pair_load_power_W','total_pair_load_power_W','world_z_mm'},forcePoints);
write_csv(fullfile(out,'material_root_points.csv'),{'source_index','interval_end_s','actual_cut_substep_s','tooth','side','material_index','s_mm','z_mm','root_height_mm','removed_height_mm','wheel_implicit_level_NOT_distance','branch','branch_gap','normal_angle_deg','target_inward_dot_vrel_mmps','wheel_outward_dot_vrel_mmps'},rootPoints);
write_csv(fullfile(out,'frame_power_and_force.csv'),{'source_index','time_s','normal_pair_load_power_W','tangent_pair_load_power_W','total_pair_load_power_W','Fx_N','Fy_N','Fz_N','Tc_Nm'},frameTotals);
zipFile=fullfile(out,'r30_contact_direction_result.zip');
Aout=struct('kind','native_replay_contact_direction_audit_v1','status','NUMERICAL_RECONSTRUCTION_PASS_PHYSICAL_INTERPRETATION_REQUIRED',...
 'outputDir',out,'resultZipFile',zipFile,'sourceIndices',selected,'selection',selection,'forceDataChanged',false,'historyRecomputed',false,...
 'forceQuadraturePointCount',size(forcePoints,1),'materialRootEventCount',size(rootPoints,1),'maxErrors',maxerr,...
 'normalPairLoadPower_W_minmax',[min(frameTotals(:,3)) max(frameTotals(:,3))],...
 'tangentPairLoadPower_W_minmax',[min(frameTotals(:,4)) max(frameTotals(:,4))],...
 'forcePointMeanAngle_FnWeighted_deg',sum(forcePoints(:,11).*forcePoints(:,13))/sum(forcePoints(:,11)),...
 'rootPointAngleMinMax_deg',[min(rootPoints(:,14)) max(rootPoints(:,14))],...
 'limitations','Force points are support-polygon centroids at interval end; root events are original material samples at substep time. Neither is a measured continuous physical contact patch. Positive reference-normal component power alone is not a proof of incorrect total force.');
save(fullfile(out,'contact_direction_audit.mat'),'Aout','-v7');write_text(fullfile(out,'contact_direction_report.json'),jsonencode(Aout));
write_text(fullfile(out,'README_zh.txt'),sprintf(['这是原生重放帧诊断，原力、材料、核心代码未改动。\n',...
 'force_quadrature_points是原条带力积分的支撑多边形重心，方向在区间末端取值。\n',...
 'material_root_points是原固定材料采样点，每次切削后的法向交点按该次子步时刻检查。\n',...
 '隐式level不是法向距离；分支1螺旋齿面、2外圆、3轴端、4根部核心。\n',...
 'normal/tangent_pair_load_power为接触载荷对双方做功，负值耗能，正值取决于参考分解的定义。\n',...
 '几何核验不要求切削点的法向相对速度等于零，因为模型允许材料去除。\n',...
 '选用帧可能是目标峰值/变号附近最近已保存帧，偏移量见selected_frames.csv。\n',...
 'PASS仅针对重建和几何恒等式，不代表模型实测精度，也不自动替换法向。\n']));
files={'selected_frames.csv','force_quadrature_points.csv','material_root_points.csv','frame_power_and_force.csv','contact_direction_audit.mat','contact_direction_report.json','README_zh.txt'};
zip(zipFile,files,out);
fprintf('\n方向诊断完成。重建最大力差%.3g N、力矩差%.3g Nmm；子步求交level残差%.3g。\n',maxerr.force_N,maxerr.moment_Nmm,maxerr.root_level);
fprintf('请上传这一份结果：\n%s\n',zipFile);
end

function validate_leaf(R,M,g)
tooth=double(M.toothIds(g));assert(isequal(R.faces,[tooth 1;tooth 2])&&isequaln(R.parameters,M.parameters)&&isequaln(R.forceCoefficients,M.coefficients));
assert(isequal(R.momentOrigin_mm,[0 0 0])&&numel(R.poses)==numel(M.time_s)&&numel(R.frameFiles)==numel(R.frameIndices));
delta=vertcat(R.poses.force_N)-reshape(M.toothForce_N(:,g,:),[],3);assert(max0(abs(delta))<1e-8);
delta=vertcat(R.poses.moment_Nmm)-reshape(M.toothMoment_Nmm(:,g,:),[],3);assert(max0(abs(delta))<1e-6);
end
function v=relative_velocity(P,A)
n=size(P,1);vw=cross(repmat(A.omegaW*A.axisW(:)',n,1),bsxfun(@minus,P,[A.center 0 0]),2);
vg=[-A.omegaG*P(:,2),A.omegaG*P(:,1),-A.feed*ones(n,1)];v=vw-vg;
end
function h=unpack(H,Q),h=permute(reshape(H,Q.n,Q.ns,Q.n,Q.nz),[2 4 1 3]);end
function x=max0(v),if isempty(v),x=0;else,x=max(v(:));end;end
function s=core_signature(core)
files=dir(fullfile(core,'*.m'));s='';for j=1:numel(files),s=[s files(j).name char(10) fileread(fullfile(core,files(j).name))];end %#ok<AGROW>
end
function out=unique_folder(base),out=base;k=0;while isfolder(out),k=k+1;out=[base '_' num2str(k)];end;end
function write_csv(path,headers,data)
fid=fopen(path,'w');assert(fid>0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',strjoin(headers,','));fmt=[repmat('%.17g,',1,size(data,2)-1) '%.17g\n'];fprintf(fid,fmt,data');
end
function write_text(path,value),fid=fopen(path,'w','n','UTF-8');assert(fid>0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',value);end
