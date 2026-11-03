function M=run_r30_multitooth_pilot(mode,user)
% RUN_R30_MULTITOOTH_PILOT 15候选齿的探索性试算入口，正式计算由本机MATLAB执行。
% 这些根目录脚本不修改matlab_recompute，不改变原双侧结果或检查点。
% 默认只筛选：
%   S=run_r30_multitooth_pilot;
% 按当前理想几何执行探索性重算：
%   M=run_r30_multitooth_pilot('run');
% 说明：当前名义外径与求交外径不同；此入口不擅自修改任何几何。
% 每个候选齿两侧从t=0独立推进；各齿分批串行求解，最后按同一时间轴
% 叠加力向量和力矩。包围界证明整个窗口不可能参与的齿直接计零贡献，
% 不为这些齿伪造材料历史。默认不保存frames，只保存曲线、状态和汇总。
if nargin<1||isempty(mode),mode='screen';end
if nargin<2,user=struct();end
assert(any(strcmp(mode,{'screen','run'})),'R30M:Mode','mode只能是screen或run。');
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'matlab_recompute'));
o=struct('parameterFile',fullfile(root,'parameters.json'),'toothIds',10:24,...
    'window_s',[14.959 15.069],'ns',80,'nz',160,'order',2,'dt_s',.0000125,...
    'substeps',2,'historyDt_s',.00005,'forceStripWidth_mm',.1,'screenStep_s',.001,...
    'loadCountThreshold_N',.01,'outputDir','','plot',true);
names=fieldnames(user);for j=1:numel(names),assert(isfield(o,names{j}),'R30M:Options','未知选项：%s',names{j});o.(names{j})=user.(names{j});end
assert(isscalar(o.loadCountThreshold_N)&&isfinite(o.loadCountThreshold_N)&&o.loadCountThreshold_N>=0);
S=r30m_preflight(o.parameterFile,o.window_s,o.toothIds,o.screenStep_s);
if strcmp(mode,'screen'),M=S;return;end
assert(S.candidateCoverageComplete,'R30M:Coverage','候选齿没有覆盖完整窗口，请扩大toothIds，不能忽略遗漏贡献。');
P=r30b_read_json(o.parameterFile);o.toothIds=o.toothIds(:)';o.window_s=o.window_s(:)';
steps=diff(o.window_s)/o.dt_s;assert(abs(steps-round(steps))<1e-7,'R30M:Window','窗口必须包含整数个力计算区间。');nt=round(steps);
fprintf('探索性试算：粗空间网格%dx%d x %d^2；力区间%.6g us；历史步长%.6g us。\n',o.ns,o.nz,o.order,o.dt_s*1e6,o.historyDt_s*1e6);
fprintf('当前几何、力系数未物理标定；12.5us的既有筛选只针对17齿局部窗口，不能视为本多齿结果已收敛。\n');
if isempty(o.outputDir)
    out=fullfile(root,['r30_multitooth_pilot_' datestr(now,'yyyymmdd_HHMMSS')]);
    assert(~isfolder(out),'R30M:Output','输出目录已存在，请指定新outputDir。');
else,out=o.outputDir;end
if ~isfolder(out),mkdir(out);end
[ok,a]=fileattrib(out);assert(ok);out=a.Name;o.outputDir=out;
planFile=fullfile(out,'pilot_plan.mat');plan=struct('parameters',P,'options',o);
if isfile(planFile)
    D=load(planFile,'plan');assert(isequaln(plan,D.plan),'R30M:Resume','续算参数不同；请使用新输出目录。');
else,save(planFile,'plan','-v7');end
write_json(fullfile(out,'screen_report.json'),S);
teeth=o.toothIds;ng=numel(teeth);time=o.window_s(1)+(1:nt)'*o.dt_s;
N=zeros(nt,3);F=N;moment=N;V=zeros(nt,1);FnScalar=V;
toothN=zeros(nt,ng,3);toothF=toothN;toothM=toothN;toothFn=zeros(nt,ng);toothVdot=toothFn;toothFnScalar=toothFn;
faceN=zeros(nt,2*ng,3);faceF=faceN;faceM=faceN;faceFn=zeros(nt,2*ng);faceVdot=faceFn;
calculated=false(1,ng);files=cell(ng,1);defaults=r30b_options;coeff=defaults.coefficients;maximumClosure=0;
for g=1:ng
    tooth=teeth(g);
    if ~ismember(tooth,S.candidateUnionTeeth)
        fprintf('齿%d：包围界已证明整个窗口不可能参与，窗口力贡献为零；不计算其材料历史。\n',tooth);
        files{g}='';continue;
    end
    fprintf('\n候选齿批次%d/%d：第%d齿，两侧独立历史。\n',g,ng,tooth);
    local=struct('faces',[tooth 1;tooth 2],'ns',o.ns,'nz',o.nz,'order',o.order,...
        'window_s',o.window_s,'dt_s',o.dt_s,'substeps',o.substeps,'historyDt_s',o.historyDt_s,...
        'forceStripWidth_mm',o.forceStripWidth_mm,'saveFrames',false,'plot',false,'visible',false,...
        'resume',true,'outputDir',fullfile(out,sprintf('tooth_%02d',tooth)));
    R=r30b_run(o.parameterFile,local);
    assert(isequal(R.faces,[tooth 1;tooth 2])&&numel(R.poses)==nt&&...
        max(abs([R.poses.t_s]'-time))<1e-10&&isequaln(R.parameters,P),...
        'R30M:Time','各齿参数、时间轴或齿面标识不一致。');
    assert(isequal(R.momentOrigin_mm,[0 0 0]),'R30M:Origin','力矩原点不一致。');
    assert(isequaln(coeff,R.forceCoefficients));
    aN=vertcat(R.poses.normalVector_N);aF=vertcat(R.poses.force_N);aM=vertcat(R.poses.moment_Nmm);aV=[R.poses.dV_mm3]';
    % 先相加向量，最后才取模，不能把各齿曲线模值直接相加当作总力。
    N=N+aN;F=F+aF;moment=moment+aM;V=V+aV;FnScalar=FnScalar+[R.poses.Fn_magnitude_sum_N]';
    toothN(:,g,:)=reshape(aN,nt,1,3);toothF(:,g,:)=reshape(aF,nt,1,3);toothM(:,g,:)=reshape(aM,nt,1,3);
    toothFn(:,g)=sqrt(sum(aN.^2,2));toothVdot(:,g)=aV/o.dt_s;toothFnScalar(:,g)=[R.poses.Fn_magnitude_sum_N]';
    for side=1:2
        q=2*(g-1)+side;ss=R.facePoses{side};faceN(:,q,:)=reshape(vertcat(ss.normalVector_N),nt,1,3);
        faceF(:,q,:)=reshape(vertcat(ss.force_N),nt,1,3);faceM(:,q,:)=reshape(vertcat(ss.moment_Nmm),nt,1,3);
        faceFn(:,q)=[ss.Fn_pure_resultant_N]';faceVdot(:,q)=[ss.Vdot_mm3ps]';
        % 独立确认所有真实去除都落在全30齿筛选保留的候选区间中。
        ds=S.timeEdges_s(2)-S.timeEdges_s(1);n=size(S.possibleFaces,1);
        lo=min(n,max(1,floor((time-o.dt_s-o.window_s(1))/ds+1e-8)+1));
        hi=min(n,max(1,ceil((time-o.window_s(1))/ds-1e-8)));lo=min(lo,hi);
        prefix=[0;cumsum(S.possibleFaces(:,2*(tooth-1)+side))];
        possible=prefix(hi+1)-prefix(lo)>0;
        assert(all(faceVdot(:,q)==0|possible),'R30M:Gate','出现筛选范围之外的去除贡献。');
    end
    maximumClosure=max(maximumClosure,max([R.validation.historyClosure_mm3(:);R.validation.windowClosure_mm3(:)]));
    calculated(g)=true;files{g}=fullfile(sprintf('tooth_%02d',tooth),'r30b_result.mat');
    clear R;
end
normalNorm=sqrt(sum(N.^2,2));fullNorm=sqrt(sum(F.^2,2));
nRemovingTeeth=sum(toothVdot>0,2);nRemovingFaces=sum(faceVdot>0,2);
nLoadedTeeth=sum(toothFnScalar>o.loadCountThreshold_N,2);
% 阈值仅用于计数显示；任何小力都没有从N/F/moment中删掉。
checkN=squeeze(sum(toothN,2));checkF=squeeze(sum(toothF,2));checkM=squeeze(sum(toothM,2));
err=max(abs([checkN(:)-N(:);checkF(:)-F(:);checkM(:)-moment(:)]));assert(err<1e-8);
M=struct('kind','exploratory_multitooth_candidate_model_not_calibrated',...
    'toothIds',teeth,'faces',zeros(2*ng,2),...
    'time_s',time,'window_s',o.window_s,'forceDt_s',o.dt_s,'geometrySubstep_s',o.dt_s/o.substeps,...
    'options',o,'parameters',P,'coefficients',coeff,'preflight',S,'calculatedHistoryForTooth',calculated,...
    'perToothResultFiles',{files},'normalVector_N',N,'force_N',F,'moment_Nmm',moment,...
    'Fn_pure_resultant_N',normalNorm,'F_resultant_N',fullNorm,'Fn_magnitude_sum_N',FnScalar,...
    'sumToothNormalNorm_N',sum(toothFn,2),'Vdot_mm3ps',V/o.dt_s,...
    'toothNormalVector_N',toothN,'toothForce_N',toothF,'toothMoment_Nmm',toothM,...
    'toothPureNormal_N',toothFn,'toothFnMagnitudeSum_N',toothFnScalar,'toothVdot_mm3ps',toothVdot,...
    'faceNormalVector_N',faceN,'faceForce_N',faceF,'faceMoment_Nmm',faceM,...
    'facePureNormal_N',faceFn,'faceVdot_mm3ps',faceVdot,...
    'nRemovingTeeth',nRemovingTeeth,'nRemovingFaces',nRemovingFaces,...
    'nLoadedTeeth_aboveThreshold',nLoadedTeeth,'loadCountThreshold_N',o.loadCountThreshold_N,...
    'numericalBookkeeping',struct('maxAggregationDifference',err,'maxMaterialClosure_mm3',maximumClosure),...
    'scope','Only above-base involute flanks with active removal; no zero-removal rubbing/contact load, no root fillet, no elastic/thermal coupling or force calibration. Original thin-normal-layer volume measure retained.');
% 所选齿面顺序为每个齿的1、2侧，明确建立而不依赖reshape推断。
M.faces=zeros(2*ng,2);for g=1:ng,M.faces(2*g-1:2*g,:)=[teeth(g) 1;teeth(g) 2];end
save(fullfile(out,'multitooth_pilot.mat'),'M','-v7');
write_summary(fullfile(out,'multitooth_summary.csv'),M);
write_matrix(fullfile(out,'per_tooth_normal_force.csv'),time,toothFn,arrayfun(@(x)sprintf('tooth%d_Fn_N',x),teeth,'UniformOutput',false));
labels=cell(1,2*ng);for f=1:2*ng,labels{f}=sprintf('tooth%d_side%d_Vdot_mm3ps',M.faces(f,:));end
write_matrix(fullfile(out,'per_face_Vdot.csv'),time,faceVdot,labels);
report=struct('status','bookkeeping_pass_exploratory_results','selectedToothIds',teeth,...
    'candidateUnionTeeth',S.candidateUnionTeeth,'candidateCoverageComplete',S.candidateCoverageComplete,...
    'window_s',o.window_s,'forceDt_s',o.dt_s,'geometrySubstep_s',o.dt_s/o.substeps,...
    'maxCandidateTeeth',S.maxCandidateTeeth,'maxRemovingTeeth',max(nRemovingTeeth),...
    'maxRemovingFaces',max(nRemovingFaces),'maxLoadedTeethAboveThreshold',max(nLoadedTeeth),...
    'countThreshold_N',o.loadCountThreshold_N,'countThresholdQuantity','sum of scalar local normal loads per tooth, before directional cancellation','thresholdChangesForceSum',false,...
    'nominalWheelDiameter_mm',S.geometry.nominalWheelOuterDiameter_mm,...
    'usedWheelDiameter_mm',S.geometry.usedWheelOuterDiameter_mm,...
    'geometryConventionAgrees',S.geometry.diameterConventionAgrees,...
    'physicalAccuracyValidated',false,'spatialConvergenceValidated',false,'multiToothTimeConvergenceValidated',false,...
    'numericalBookkeeping',M.numericalBookkeeping);
write_json(fullfile(out,'pilot_report.json'),report);
upload={'multitooth_pilot.mat','multitooth_summary.csv','per_tooth_normal_force.csv','per_face_Vdot.csv','pilot_report.json','screen_report.json'};
if o.plot
    try
        fig=figure('Visible','off','Color','w');guard=onCleanup(@()close(fig));
        subplot(5,1,1);imagesc(time,1:ng,toothFn');axis xy;yticks(1:ng);yticklabels(arrayfun(@num2str,teeth,'UniformOutput',false));ylabel('Tooth');colorbar;title('Exploratory per-tooth pure normal resultant / N');
        subplot(5,1,2);stairs(time,nRemovingTeeth);hold on;stairs(time,nRemovingFaces);stairs(time,nLoadedTeeth);grid on;legend('Removing teeth','Removing faces','Teeth above display threshold');ylabel('Count');
        subplot(5,1,3);plot(time,[normalNorm sum(toothFn,2)]);grid on;legend('Vector-sum norm','Sum of tooth norms');ylabel('Pure normal / N');
        subplot(5,1,4);plot(time,F);grid on;legend('Fx','Fy','Fz');ylabel('Full force / N');
        subplot(5,1,5);plot(time,moment(:,3)/1000);grid on;ylabel('C moment / Nm');xlabel('Time / s');
        saveas(fig,fullfile(out,'multitooth_pilot.png'));upload{end+1}='multitooth_pilot.png';clear guard;
    catch e,write_text(fullfile(out,'graphics_error.txt'),e.message);fprintf('数值结果已保存，绘图未完成：%s\n',e.message);end
end
zipFile=fullfile(out,'multitooth_pilot_summary.zip');zip(zipFile,upload,out);
fprintf('探索性试算完成：最大实际去除齿数%d，最大去除齿面数%d。\n',max(nRemovingTeeth),max(nRemovingFaces));
fprintf('请提供汇总包：%s\n',zipFile);
end

function write_summary(path,M)
values=[M.time_s M.force_N M.normalVector_N M.Fn_pure_resultant_N M.sumToothNormalNorm_N ...
    M.Fn_magnitude_sum_N M.moment_Nmm(:,3)/1000 M.Vdot_mm3ps M.nRemovingTeeth M.nRemovingFaces M.nLoadedTeeth_aboveThreshold];
labels={'Fx_N','Fy_N','Fz_N','Nx_N','Ny_N','Nz_N','Fn_vector_sum_N','sum_tooth_Fn_norm_N',...
    'sum_strip_Fn_N','Tc_Nm','Vdot_mm3ps','nRemovingTeeth','nRemovingFaces','nLoadedTeeth_aboveThreshold'};
write_matrix(path,M.time_s,values(:,2:end),labels);
end
function write_matrix(path,time,values,labels)
fid=fopen(path,'w');assert(fid>0);guard=onCleanup(@()fclose(fid));fprintf(fid,'t_s');for j=1:numel(labels),fprintf(fid,',%s',labels{j});end;fprintf(fid,'\n');
for k=1:numel(time),v=[time(k) values(k,:)];fprintf(fid,'%.16g',v(1));fprintf(fid,',%.16g',v(2:end));fprintf(fid,'\n');end
end
function write_json(path,value),write_text(path,jsonencode(value));end
function write_text(path,text)
fid=fopen(path,'w','n','UTF-8');assert(fid>0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',text);
end
