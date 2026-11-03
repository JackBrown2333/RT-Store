function S=r30m_preflight(parameterFile,window_s,toothIds,screenStep_s)
% R30M_PREFLIGHT 全30齿的保守包围界筛选，不计算材料去除或磨削力。
% 放在run_r30_dual.m同级： S=r30m_preflight;
% possible表示可能碰到包围砂轮，不等于真实接触、去除或承载。
if nargin<1||isempty(parameterFile),parameterFile=fullfile(fileparts(mfilename('fullpath')),'parameters.json');end
if nargin<2||isempty(window_s),window_s=[14.959 15.069];end
if nargin<3||isempty(toothIds),toothIds=10:24;end
if nargin<4||isempty(screenStep_s),screenStep_s=.001;end
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'matlab_recompute'));
P=r30b_read_json(parameterFile);window_s=window_s(:)';toothIds=toothIds(:)';
assert(numel(window_s)==2&&all(isfinite(window_s))&&window_s(1)>=0&&window_s(2)>window_s(1));
assert(all(isfinite(toothIds))&&all(ismember(toothIds,1:P.p.gear.teeth))&&numel(unique(toothIds))==numel(toothIds));
assert(isscalar(screenStep_s)&&isfinite(screenStep_s)&&screenStep_s>0);
n=P.p.gear.teeth;faces=zeros(2*n,2);
for t=1:n,faces(2*t-1:2*t,:)=[t 1;t 2];end
edges=linspace(window_s(1),window_s(2),ceil(diff(window_s)/screenStep_s)+1);
possible=false(numel(edges)-1,2*n);
for k=1:size(possible,1)
    for f=1:size(faces,1),possible(k,f)=r30b_gate(P,faces(f,:),edges(k),edges(k+1),1);end
end
toothPossible=false(size(possible,1),n);
for t=1:n,toothPossible(:,t)=any(possible(:,2*t-1:2*t),2);end
unionTeeth=find(any(toothPossible,1));missing=setdiff(unionTeeth,toothIds);
diameterUsed=2*P.installation.tipNew;nominal=P.p.calc.wheel.outerDiameter_mm;
S=struct('kind','conservative_candidate_screen_not_actual_engagement','faces',faces,...
    'window_s',window_s,'timeEdges_s',edges,'possibleFaces',possible,'possibleTeeth',toothPossible,...
    'selectedToothIds',toothIds,'candidateUnionTeeth',unionTeeth,'missingCandidateTeeth',missing,...
    'candidateCoverageComplete',isempty(missing),'maxCandidateTeeth',max(sum(toothPossible,2)),...
    'maxCandidateFaces',max(sum(possible,2)),'normalLayerBound_mm',[-1 P.cfg.stock_um/1000],...
    'geometry',struct('nominalWheelOuterDiameter_mm',nominal,'usedWheelOuterRadius_mm',P.installation.tipNew,...
    'usedWheelOuterDiameter_mm',diameterUsed,'diameterConventionAgrees',abs(diameterUsed-nominal)<1e-8,...
    'profilePhysicallyValidated',false),...
    'note','覆盖整个筛选小时间区间的保守几何上界；未重建材料历史，不能将候选数量当作实际接触或承载齿数。');
fprintf('共同窗口：%.6f..%.6f s\n',window_s);
fprintf('全窗口候选齿：');fprintf('%d ',unionTeeth);fprintf('\n');
fprintf('任一筛选区间候选上界：%d齿、%d齿面。\n',S.maxCandidateTeeth,S.maxCandidateFaces);
fprintf('名义外径 %.6f mm；当前求交外径 %.6f mm。\n',nominal,diameterUsed);
if ~S.geometry.diameterConventionAgrees
    fprintf('注意：这是两种几何口径；正式物理预测前需确认实际砂轮型面和外径。\n');
end
if ~isempty(missing),fprintf('所选齿遗漏候选：');fprintf('%d ',missing);fprintf('\n');end
end
