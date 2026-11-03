function zipFiles=r30m_export_multitooth_animation(resultDir,user)
% R30M_EXPORT_MULTITOOTH_ANIMATION 从多齿窗口初态重放并导出动画数据。
% 放在run_r30_multitooth_pilot.m同级；不改matlab_recompute或原结果。
%   zipFiles=r30m_export_multitooth_animation;
% 默认读取r30_multitooth_pilot_15_v2_160x320，输出到其animation_export子目录。
% 力曲线直接复制；保留帧的几何与条带复算后逐项核对。每个时间子步
% 都会推进材料，抽帧之间的材料变化累计；不重新计算窗口前的历史。
if nargin<1||isempty(resultDir)
    resultDir=fullfile(fileparts(mfilename('fullpath')),'r30_multitooth_pilot_15_v2_160x320');
end
if nargin<2,user=struct();end
root=fileparts(mfilename('fullpath'));core=fullfile(root,'matlab_recompute');addpath(core);
assert(isfile(fullfile(root,'r30b_export_animation.m')),'R30M:Dependency','请将本包的r30b_export_animation.m也放在工程根目录。');
[ok,a]=fileattrib(resultDir);assert(ok&&a.directory,'R30M:Source','找不到多齿结果目录。');sourceDir=a.Name;
sourceMat=fullfile(sourceDir,'multitooth_pilot.mat');assert(isfile(sourceMat),'R30M:Source','缺少multitooth_pilot.mat。');D=load(sourceMat,'M');M=D.M;clear D;
o=struct('frameStride',40,'profileSamples',96,'axialSamples',128,'framesPerChunk',16,...
    'targetPartMiB',15,'keepCriticalFrames',true,'outputDir',fullfile(sourceDir,'animation_export'),...
    'progressEvery_s',5,'checkpointEvery_s',30);
names=fieldnames(user);for j=1:numel(names),assert(isfield(o,names{j}),'R30M:Options','未知选项：%s',names{j});o.(names{j})=user.(names{j});end
v=[o.frameStride o.profileSamples o.axialSamples o.framesPerChunk];assert(all(isfinite(v))&&all(v>=1)&&all(mod(v,1)==0));
assert(isscalar(o.targetPartMiB)&&isfinite(o.targetPartMiB)&&o.targetPartMiB>=.1);
assert(o.progressEvery_s>0&&o.checkpointEvery_s>0&&isscalar(o.keepCriticalFrames));
nt=numel(M.time_s);ng=numel(M.toothIds);assert(nt>0&&size(M.faces,1)==2*ng);
keep=select_frames(M,o);coreSignature=source_signature(core);
for pair={{'toothForce_N','force_N'},{'toothNormalVector_N','normalVector_N'},{'toothMoment_Nmm','moment_Nmm'}}
    names=pair{1};totalCheck=reshape(sum(M.(names{1}),2),nt,3)-M.(names{2});
    assert(max(abs(totalCheck(:)))<1e-7,'R30M:Summary','逐齿力或力矩汇总不符。');
end
assert(max(abs(sum(M.toothVdot_mm3ps,2)-M.Vdot_mm3ps))<1e-8,'R30M:Summary','逐齿去除率汇总不符。');
assert(isequal(sum(M.toothVdot_mm3ps>0,2),double(M.nRemovingTeeth))&&...
    isequal(sum(M.faceVdot_mm3ps>0,2),double(M.nRemovingFaces)),'R30M:Summary','去除数量与体积率不符。');
if ~isfolder(o.outputDir),mkdir(o.outputDir);end
[ok,a]=fileattrib(o.outputDir);assert(ok);out=a.Name;o.outputDir=out;
plan=struct('schema','r30_multitooth_replay_v1','sourceDir',sourceDir,'M',M,'options',o,...
    'selectedIndices',keep,'coreSignature',coreSignature,...
    'replaySource',fileread([mfilename('fullpath') '.m']),'compactExporterSource',fileread(fullfile(root,'r30b_export_animation.m')));
planFile=fullfile(out,'export_plan.mat');
if isfile(planFile)
    D=load(planFile,'plan');assert(isequaln(D.plan,plan),'R30M:Resume','源结果、代码或导出设置已改变，请指定新的outputDir。');clear D;
else
    listing=dir(out);assert(numel(listing)<=2,'R30M:Output','输出目录非空且没有本工具计划文件，请指定新目录。');
    save(planFile,'plan','-v7');
end
clear plan;
fprintf('动画准备：%d候选齿，%d个原力区间，保留%d个共同几何帧。\n',ng,nt,numel(keep));
fprintf('从保存的窗口初态重放 %.6f..%.6f s；不重跑窗口前历史，曲线沿用原结果。\n',M.window_s);
entries={};checks={};
for g=1:ng
    tooth=double(M.toothIds(g));
    if ~M.calculatedHistoryForTooth(g)
        assert(all(M.toothVdot_mm3ps(:,g)==0)&&all(reshape(M.toothForce_N(:,g,:),[],1)==0),'R30M:Skipped','未算历史的齿存在非零贡献。');
        continue;
    end
    relative=strrep(M.perToothResultFiles{g},'\',filesep);leaf=fullfile(sourceDir,relative);leafDir=fileparts(leaf);
    assert(isfile(leaf),'R30M:Missing','缺少原始逐齿结果：%s。请使用本机完整结果目录。',leaf);
    D=load(leaf,'R30B');R=D.R30B;clear D;validate_leaf(R,M,g);
    cp=fullfile(leafDir,'r30b_checkpoint.mat');assert(isfile(cp),'R30M:Missing','缺少逐齿检查点：%s',cp);
    D=load(cp,'S');assert(strcmp(D.S.phase,'done')&&strcmp(D.S.meta.codeSignature,coreSignature),...
        'R30M:Core','原逐齿计算未完成或matlab_recompute代码已改变，请恢复计算时使用的核心代码。');clear D;
    D=load(fullfile(leafDir,R.initialStateFile),'initialState');initial=D.initialState;clear D;
    D=load(fullfile(leafDir,R.finalStateFile),'finalState');final=D.finalState;clear D;
    working=fullfile(out,'replay',sprintf('tooth_%02d',tooth));if ~isfolder(working),mkdir(working);end
    verify_source_copies(working,R,initial,final);
    fprintf('\n重放第%d齿的两个齿面（批次%d/%d）。\n',tooth,g,ng);
    c=replay_leaf(R,initial,final,working,keep,o);checks{end+1}=c; %#ok<AGROW>
    doneFile=fullfile(working,'compact_result.mat');reuse=false;
    if isfile(doneFile)
        D=load(doneFile,'entry');entry=D.entry;clear D;
        reuse=isfile(fullfile(entry.directory,'export_report.json'))&&all(cellfun(@(x)isfile(fullfile(entry.directory,x)),entry.files));
    end
    if ~reuse
        compact=unique_folder(fullfile(out,'compact',sprintf('tooth_%02d',tooth)));
        eo=struct('frameStride',1,'profileSamples',o.profileSamples,'axialSamples',o.axialSamples,...
            'framesPerChunk',o.framesPerChunk,'targetPartMiB',o.targetPartMiB,'keepCriticalFrames',false,'outputDir',compact);
        fprintf('  精简第%d齿的显示帧与接触边界。\n',tooth);
        compactLog=evalc('r30b_export_animation(working,eo);');
        write_text(fullfile(working,'compact_export_log.txt'),compactLog);
        fprintf('  第%d齿精简回放核验通过。\n',tooth);
        manifest=jsondecode(fileread(fullfile(compact,'animation_manifest.json')));
        chunks=to_cellstr(manifest.chunks);entry=struct('tooth',tooth,'faces',R.faces,'directory',compact,...
            'files',{{'animation_header.mat','animation_curves.mat','export_report.json',chunks{:}}},'chunks',{chunks});
        save(doneFile,'entry','-v7');
    end
    entries{end+1}=entry; %#ok<AGROW>
    clear R initial final;
end
assert(~isempty(entries),'R30M:Empty','没有已计算材料历史的齿。');
package=fullfile(out,'package');if ~isfolder(package),mkdir(package);end
copyfile(sourceMat,fullfile(package,'multitooth_pilot.mat'));
allFiles={'multitooth_pilot.mat'};bundles=cell(size(entries));
for g=1:numel(entries)
    e=entries{g};prefix=sprintf('tooth_%02d',e.tooth);dest=fullfile(package,prefix);if ~isfolder(dest),mkdir(dest);end
    for j=1:numel(e.files),copyfile(fullfile(e.directory,e.files{j}),fullfile(dest,e.files{j}));allFiles{end+1}=fullfile(prefix,e.files{j});end %#ok<AGROW>
    bundles{g}=struct('tooth',e.tooth,'faces',e.faces,'prefix',prefix,'chunks',{e.chunks});
end
report=struct('status','PASS','kind','verified_window_replay_and_display_export_not_physical_validation',...
    'sourceHistoryRecomputed',false,'sourceForceCurvesChanged',false,'selectedFrameStripsRecomputedAndChecked',true,...
    'sourceGrid',[M.options.ns M.options.nz],'sourceForceDt_s',M.forceDt_s,...
    'sourceGeometrySubstep_s',M.geometrySubstep_s,'sourceIntervals',nt,'exportedCommonFrames',numel(keep),...
    'materialTeeth',cellfun(@(e)e.tooth,entries),'perToothChecks',{checks},...
    'physicalAccuracyValidated',false,'geometryDiameter_mm',2*M.parameters.installation.tipNew);
manifest=struct('schemaVersion','r30_multitooth_animation_v1','sourceSummary','multitooth_pilot.mat',...
    'selectedToothIds',M.toothIds,'faces',M.faces,'window_s',M.window_s,...
    'selectedSourceIntervalIndices',keep,'forceDataChanged',false,'forceOutputDt_s',M.forceDt_s,...
    'materialBundles',{bundles},'historyNotComputedForTeeth',M.toothIds(~M.calculatedHistoryForTooth),...
    'scope','Computed above-base material flanks only; other selected teeth are context, not computed remaining material.',...
    'parts',{{}},'partFiles',{{}});
write_json(fullfile(package,'multitooth_animation_report.json'),report);
write_text(fullfile(package,'README_zh.txt'),sprintf(['本包是160x320等原空间网格的多齿窗口重放导出。\n',...
    '根目录multitooth_pilot.mat的M为原曲线，未缩放或拼接。\n',...
    '以multitooth_animation_manifest.json组织各tooth_XX紧凑帧；MATLAB索引从1开始。\n',...
    '每个原时间子步均推进材料；仅抽取显示帧，帧间材料更新已累计。\n',...
    '力为原时间区间平均量，帧接触轮廓属于该原区间，不代表抽帧间所有去除的并集。\n',...
    '显示表面为选取的原材料点高度；点间插值仅用于显示，接触边界为活动支撑矩形并集。\n',...
    'world尺寸锚点与保留帧条带几何重放并核对，非精确连续实体接触曲线。\n',...
    '未计算历史的齿仅作几何背景，不得冒充已计算的剩余材料。\n',...
    '当前理想砂轮边界外径与绝对力的物理标定限制保持不变。\n']));
shared={'multitooth_animation_manifest.json','multitooth_animation_report.json','README_zh.txt'};
groups=group_files(package,allFiles,o.targetPartMiB);names=cell(numel(groups),1);
for j=1:numel(groups),names{j}=sprintf('r30_multitooth_animation_%02d.zip',j);end
manifest.parts=names;manifest.partFiles=groups;write_json(fullfile(package,shared{1}),manifest);
zipFiles=cell(numel(groups),1);
for j=1:numel(groups)
    zipFiles{j}=fullfile(out,names{j});zip(zipFiles{j},[shared groups{j}],package);
    x=dir(zipFiles{j});assert(x.bytes<=o.targetPartMiB*1024^2,'R30M:PartSize','分包过大，请增加targetPartMiB并使用新outputDir。');
end
fprintf('\n多齿动画导出核验：PASS。请上传下列全部ZIP；不需要上传本机replay缓存或原检查点。\n');
for j=1:numel(zipFiles),x=dir(zipFiles{j});fprintf('%s  (%.2f MiB)\n',zipFiles{j},x.bytes/1024^2);end
end

function c=replay_leaf(R,initial,final,out,keep,eo)
o=R.options;o.saveFrames=true;o.frameStride=eo.frameStride;o.outputDir=out;
nf=size(R.faces,1);nt=numel(R.poses);assert(nf==2);
assert(abs(initial.t_s-R.window_s(1))<1e-10&&abs(final.t_s-R.window_s(2))<1e-10);
cp=fullfile(out,'animation_replay_checkpoint.mat');
if isfile(cp)
    D=load(cp,'S');S=D.S;clear D;
else
    S=struct('phase','window','next',1,'Q',{{}},'lastFine',{{}},'frameFiles',{{}},'frameIndices',[],...
        'lastFrameTime',R.window_s(1),'maxStepVolumeDifference_mm3',0,'maxSubstepBalance_mm3',0,...
        'maxSelectedForceDifference_N',0,'maxSelectedMomentDifference_Nmm',0,'maxSelectedAreaDifference_mm2',0,...
        'maxFinalHeightDifference_mm',0,'maxPassDepthDifference_mm',0,'selectedForceChecks',0);
    for f=1:nf
        q=initial.faces{f};validate_state_label(q,R.faces(f,:));validate_state_label(final.faces{f},R.faces(f,:));
        Q=r30b_init(R.parameters,o,R.faces(f,:));
        assert(isequal(size(q.height_mm),[Q.ns*Q.n Q.nz*Q.n])&&all(isfinite(q.height_mm(:))));
        assert(abs(q.stock_nominal_mm-Q.stock)<1e-14);
        assert(max(abs(q.support_s_mm(:)-Q.support_s(:)))<1e-13&&max(abs(q.support_z_mm(:)-Q.support_z(:)))<1e-13);
        ss=reshape(Q.s',[],1);zz=reshape(Q.z',[],1);
        assert(max(abs(ss-q.sample_s_mm(:)))<1e-13&&max(abs(zz-q.sample_z_mm(:)))<1e-13);
        Q.h=unpack_heights(q.height_mm,Q);Q.passStart=unpack_heights(q.passStartHeight_mm,Q);Q.passId=q.passId;
        % 原历史最后一次处理的正是窗口起点，其inPass与该点包围界一致。
        Q.inPass=r30b_gate(R.parameters,R.faces(f,:),R.window_s(1),R.window_s(1),o.normalDepthLimit_mm);
        Q.historyVolume_mm3=R.historyVolume_mm3(f);S.Q{f}=Q;S.lastFine{f}=q.height_mm;
    end
    r30b_atomic(cp,'S',S);
end
if ~isfolder(fullfile(out,'frames')),mkdir(fullfile(out,'frames'));end
clock=tic;checkpoint=tic;
for k=S.next:nt
    ta=R.window_s(1)+(k-1)*o.dt_s;t=R.window_s(1)+k*o.dt_s;capture=ismember(k,keep);
    frame=struct('engineVersion',R.engineVersion,'index',k,'t_start_s',ta,'t_s',t,'dt_s',o.dt_s,...
        'toolSampleTimes_s',ta+(1:o.substeps)*(t-ta)/o.substeps,'faces',{{}});
    for f=1:nf
        Q=S.Q{f};before=r30b_heights(Q);subVolume=0;subRecords=cell(o.substeps,1);
        for sub=1:o.substeps
            ts=ta+(t-ta)*sub/o.substeps;
            if r30b_gate(R.parameters,R.faces(f,:),ts,ts,o.normalDepthLimit_mm)
                if ~Q.inPass,Q.passStart=Q.h;Q.passId=Q.passId+1;Q.inPass=true;end
                if capture,[Q,v,~,change]=r30b_step(Q,ts,o);subRecords{sub}=change;
                else,[Q,v]=r30b_step(Q,ts,o);end
                subVolume=subVolume+v;
            else
                Q.inPass=false;
                if capture,subRecords{sub}=struct('t_s',ts,'index',zeros(0,1),'before_mm',zeros(0,1),'after_mm',zeros(0,1));end
            end
        end
        after=r30b_heights(Q);dh=before-after;volume=sum(Q.weightsFine(:).*dh(:));
        target=R.facePoses{f}(k);assert(abs(target.t_s-t)<1e-10);
        dv=abs(volume-target.dV_mm3);S.maxStepVolumeDifference_mm3=max(S.maxStepVolumeDifference_mm3,dv);
        S.maxSubstepBalance_mm3=max(S.maxSubstepBalance_mm3,abs(volume-subVolume));
        assert(dv<1e-9&&abs(volume-subVolume)<1e-9,'R30M:ReplayVolume','第%d齿第%d侧第%d区间体积重放不符，已停止导出。',R.faces(f,:),k);
        if capture
            temp=Q;temp.h=Q.passStart;pass=r30b_heights(temp);
            rows=r30b_force(r30b_strips(Q,before,after,pass,t,o.dt_s,o,true),o.coefficients);
            area=sum(Q.weightsFine(dh>0));s=r30b_sum(rows,o.referenceOutward,t,o.dt_s,volume,area);
            ef=max(abs([s.force_N-target.force_N s.normalVector_N-target.normalVector_N]));em=max(abs(s.moment_Nmm-target.moment_Nmm));ea=abs(area-target.area_mm2);
            S.maxSelectedForceDifference_N=max(S.maxSelectedForceDifference_N,ef);S.maxSelectedMomentDifference_Nmm=max(S.maxSelectedMomentDifference_Nmm,em);S.maxSelectedAreaDifference_mm2=max(S.maxSelectedAreaDifference_mm2,ea);
            ep=abs(s.ae_pass_peak_mm-target.ae_pass_peak_mm);S.maxPassDepthDifference_mm=max(S.maxPassDepthDifference_mm,ep);
            assert(ef<1e-7&&em<1e-5&&ea<1e-8&&ep<1e-9,'R30M:ReplayForce','保留帧的力、力矩、面积或扫掠深度与原结果不符。');
            S.selectedForceChecks=S.selectedForceChecks+1;
            ids=find(dh>0);[i,j]=ind2sub(size(after),ids);changed=find(after<S.lastFine{f});
            removal=struct('index',ids,'before_mm',before(ids),'after_mm',after(ids),'passStart_mm',pass(ids));
            frame.faces{f}=struct('tooth_id',R.faces(f,1),'flank_side',R.faces(f,2),'passId',Q.passId,...
                'stateUpdates',struct('since_t_s',S.lastFrameTime,'index',changed,'height_mm',after(changed)),...
                'removal',removal,'substepUpdates',{subRecords},'strips',rows,'summary',target);
            S.lastFine{f}=after;
        end
        S.Q{f}=Q;
    end
    if capture
        frame.combined=R.poses(k);name=sprintf('frames/frame_%06d.mat',k);r30b_atomic(fullfile(out,name),'frame',frame);
        S.frameFiles{end+1}=name;S.frameIndices(end+1)=k;S.lastFrameTime=t;
    end
    S.next=k+1;
    if toc(clock)>eo.progressEvery_s,fprintf('  窗口重放 %d/%d，已保存%d帧。\n',k,nt,numel(S.frameIndices));clock=tic;end
    if toc(checkpoint)>eo.checkpointEvery_s,r30b_atomic(cp,'S',S);checkpoint=tic;end
end
for f=1:nf
    H=r30b_heights(S.Q{f});expected=final.faces{f}.height_mm;assert(isequal(size(H),size(expected)));
    e=max(abs(H(:)-expected(:)));S.maxFinalHeightDifference_mm=max(S.maxFinalHeightDifference_mm,e);
end
assert(S.maxFinalHeightDifference_mm<1e-12,'R30M:ReplayFinal','末态材料高度与原结果不一致。');
S.phase='done';r30b_atomic(cp,'S',S);
R.frameFiles=S.frameFiles;R.frameIndices=S.frameIndices;R.options=o;
R.animationReplayOnly=true;R.animationHistoryRecomputed=false;
R30B=R;save(fullfile(out,'r30b_result.mat'),'R30B','-v7');
c=struct('tooth',R.faces(1,1),'status','PASS','sourceHistoryFromZero',R.historyFromZero,...
    'historyRecomputed',false,'forceCurvesChanged',false,'sourceIntervals',nt,'selectedFrames',numel(S.frameIndices),...
    'maxStepVolumeDifference_mm3',S.maxStepVolumeDifference_mm3,'maxSubstepBalance_mm3',S.maxSubstepBalance_mm3,...
    'maxSelectedForceDifference_N',S.maxSelectedForceDifference_N,'maxSelectedMomentDifference_Nmm',S.maxSelectedMomentDifference_Nmm,...
    'maxSelectedAreaDifference_mm2',S.maxSelectedAreaDifference_mm2,'maxPassDepthDifference_mm',S.maxPassDepthDifference_mm,...
    'maxFinalHeightDifference_mm',S.maxFinalHeightDifference_mm,'selectedForceChecks',S.selectedForceChecks);
write_json(fullfile(out,'animation_replay_report.json'),c);
end

function h=unpack_heights(H,Q),h=permute(reshape(H,Q.n,Q.ns,Q.n,Q.nz),[2 4 1 3]);end
function validate_state_label(q,face),assert(q.tooth_id==face(1)&&q.flank_side==face(2),'R30M:State','材料状态齿面标签不符。');end
function validate_leaf(R,M,g)
tooth=double(M.toothIds(g));assert(strcmp(R.engineVersion,'r30_dual_from_zero_v2')&&R.historyFromZero&&~R.oldResultsUsed);
assert(isequal(R.faces,[tooth 1;tooth 2])&&isequaln(R.parameters,M.parameters)&&isequaln(R.forceCoefficients,M.coefficients));
assert(isequal(R.momentOrigin_mm,[0 0 0])&&numel(R.poses)==numel(M.time_s));
fields={'ns','nz','order','dt_s','substeps','historyDt_s','forceStripWidth_mm','window_s'};
for j=1:numel(fields),assert(isequaln(R.options.(fields{j}),M.options.(fields{j})),'R30M:Grid','逐齿源设置与多齿汇总不同。');end
assert(max(abs([R.poses.t_s]'-M.time_s))<1e-10);
assert(max(abs(reshape(vertcat(R.poses.force_N)-reshape(M.toothForce_N(:,g,:),[],3),[],1)))<1e-8);
assert(max(abs(reshape(vertcat(R.poses.normalVector_N)-reshape(M.toothNormalVector_N(:,g,:),[],3),[],1)))<1e-8);
assert(max(abs(reshape(vertcat(R.poses.moment_Nmm)-reshape(M.toothMoment_Nmm(:,g,:),[],3),[],1)))<1e-6);
assert(max(abs([R.poses.Vdot_mm3ps]'-M.toothVdot_mm3ps(:,g)))<1e-8);
for side=1:2
    j=2*g-2+side;s=R.facePoses{side};
    assert(max(abs(reshape(vertcat(s.force_N)-reshape(M.faceForce_N(:,j,:),[],3),[],1)))<1e-8);
    assert(max(abs(reshape(vertcat(s.normalVector_N)-reshape(M.faceNormalVector_N(:,j,:),[],3),[],1)))<1e-8);
    assert(max(abs(reshape(vertcat(s.moment_Nmm)-reshape(M.faceMoment_Nmm(:,j,:),[],3),[],1)))<1e-6);
    assert(max(abs([s.Vdot_mm3ps]'-M.faceVdot_mm3ps(:,j)))<1e-8);
end
end
function verify_source_copies(out,R,initialState,finalState)
paths={'source_result.mat','initial_state.mat','final_state.mat'};names={'R30B','initialState','finalState'};values={R,initialState,finalState};
for j=1:3
    path=fullfile(out,paths{j});
    if isfile(path),D=load(path,names{j});assert(isequaln(D.(names{j}),values{j}),'R30M:Resume','源结果或材料状态改变，请换新的导出目录。');
    else,D=struct();D.(names{j})=values{j};save(path,'-struct','D','-v7');end
end
end
function indices=select_frames(M,o)
nt=numel(M.time_s);indices=unique([1:o.frameStride:nt nt]);if ~o.keepCriticalFrames,return;end
wanted=[1 nt find(diff(M.nRemovingTeeth)~=0)'+1 find(diff(M.nRemovingFaces)~=0)'+1];
for values={M.Fn_pure_resultant_N,M.F_resultant_N,abs(M.moment_Nmm(:,3))}
    [~,k]=max(values{1});wanted(end+1)=k; %#ok<AGROW>
end
for f=1:size(M.faceVdot_mm3ps,2)
    active=find(M.faceVdot_mm3ps(:,f)>0);if isempty(active),continue;end
    [~,k]=max(M.facePureNormal_N(:,f));wanted=[wanted active(1) active(end) k]; %#ok<AGROW>
    near=active(1):min(nt,active(1)+ceil(.003/M.forceDt_s));[~,j]=max(M.facePureNormal_N(near,f));wanted(end+1)=near(j); %#ok<AGROW>
end
indices=unique([indices wanted wanted-1 wanted+1]);indices=indices(indices>=1&indices<=nt);
end
function s=source_signature(core)
files=dir(fullfile(core,'*.m'));s='';for j=1:numel(files),s=[s files(j).name char(10) fileread(fullfile(core,files(j).name))];end %#ok<AGROW>
end
function out=unique_folder(base)
out=base;k=0;while isfolder(out),k=k+1;out=[base '_attempt' num2str(k)];end
end
function c=to_cellstr(x)
if iscell(x),c=x(:)';elseif ischar(x),c={x};else,c=cellstr(x);c=c(:)';end
end
function groups=group_files(root,files,target)
budget=target*1024^2*.88;groups={};batch={};bytes=0;
for j=1:numel(files)
    d=dir(fullfile(root,files{j}));assert(numel(d)==1&&d.bytes<budget,'R30M:PartSize','单个数据文件过大，请增加targetPartMiB。');
    if bytes+d.bytes>budget&&~isempty(batch),groups{end+1}=batch;batch={};bytes=0;end %#ok<AGROW>
    batch{end+1}=files{j};bytes=bytes+d.bytes; %#ok<AGROW>
end
if ~isempty(batch),groups{end+1}=batch;end
end
function write_json(path,value),write_text(path,jsonencode(value));end
function write_text(path,value)
fid=fopen(path,'w','n','UTF-8');assert(fid>0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',value);
end
