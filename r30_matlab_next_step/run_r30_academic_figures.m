function P=run_r30_academic_figures(resultDir,user)
% RUN_R30_ACADEMIC_FIGURES Saved multitooth forces -> editable paper figures.
% Put this file and r30f_spectrum.m in the project root, or run the supplied
% package directly (example_data contains the unchanged latest M snapshot).
%   P=run_r30_academic_figures;
%   P=run_r30_academic_figures('r30_multitooth_pilot_15_v2_160x320');
%   P=run_r30_academic_figures([],struct('singleTooth',17,'pairTeeth',[17 18]));
% Requires MATLAB R2020a+ for exportgraphics; no Signal Processing Toolbox.
% No force scaling, smoothing, resampling or force-model changes are made.
if nargin<1,resultDir=[];end;if nargin<2,user=struct();end
root=fileparts(mfilename('fullpath'));
projectRoot=root;if ~isfolder(fullfile(root,'matlab_recompute'))&&isfolder(fullfile(fileparts(root),'matlab_recompute')),projectRoot=fileparts(root);end
o=struct('singleTooth',17,'pairTeeth',[17 18],'fftMaxHz',500,'outputDir','',...
 'visible',true,'makeFigures',true,'allActiveTeeth',false,'fontName','Times New Roman',...
 'chineseFontName','SimSun','lineWidth',.65,'fontSize',10.5,'pngDpi',400);
names=fieldnames(user);for j=1:numel(names),assert(isfield(o,names{j}),'R30F:Option','Unknown option: %s',names{j});o.(names{j})=user.(names{j});end
assert(isscalar(o.fftMaxHz)&&isfinite(o.fftMaxHz)&&o.fftMaxHz>0);
assert(numel(o.pairTeeth)==2&&numel(unique(o.pairTeeth))==2,'R30F:Pair','Choose two distinct teeth.');
if isempty(resultDir)
 candidates={fullfile(projectRoot,'r30_multitooth_pilot_15_v2_160x320'),fullfile(root,'example_data')};
 for j=1:numel(candidates),if isfile(fullfile(candidates{j},'multitooth_pilot.mat')),resultDir=candidates{j};break;end;end
end
assert(~isempty(resultDir),'R30F:Input','Specify the folder containing multitooth_pilot.mat.');
if isfile(resultDir),sourceFile=resultDir;else,sourceFile=fullfile(resultDir,'multitooth_pilot.mat');end
assert(isfile(sourceFile),'R30F:Input','Missing multitooth_pilot.mat: %s',sourceFile);
D=load(sourceFile,'M');assert(isfield(D,'M'));M=D.M;clear D;
t=double(M.time_s(:));dt=double(M.forceDt_s);n=numel(t);ids=double(M.toothIds(:)');
assert(n>3&&all(diff(t)>0)&&max(abs(diff(t)-dt))<1e-10,'R30F:Sampling','FFT requires the original uniform time grid.');
assert(abs(n*dt-diff(M.window_s))<1e-9,'R30F:Window','Record duration and saved window disagree.');
assert(size(M.force_N,1)==n&&size(M.force_N,2)==3&&all(isfinite(M.force_N(:))));
toothF=double(M.toothForce_N);toothN=double(M.toothNormalVector_N);
active=find(any(M.toothVdot_mm3ps>0,1));
g=find(ids==o.singleTooth);pair=arrayfun(@(x)find(ids==x,1),o.pairTeeth,'UniformOutput',false);
assert(isscalar(g)&&all(cellfun(@numel,pair)==1),'R30F:Tooth','Requested teeth must be present in the saved result.');
pair=cell2mat(pair);assert(all(ismember([g pair],active)),'R30F:Inactive','Choose teeth that remove material in this window.');
F1=reshape(toothF(:,pair(1),:),n,3);F2=reshape(toothF(:,pair(2),:),n,3);Fpair=F1+F2;
Fsingle=reshape(toothF(:,g,:),n,3);Nsingle=reshape(toothN(:,g,:),n,3);
norm1=vecnorm(F1,2,2);norm2=vecnorm(F2,2,2);normPair=vecnorm(Fpair,2,2);upper=norm1+norm2;
eta=NaN(n,1);loaded=upper>1e-12;eta(loaded)=1-normPair(loaded)./upper(loaded);
assert(all(eta(loaded)>=-1e-12&eta(loaded)<=1+1e-12),'R30F:Triangle','Vector-sum bound failed.');
delta=reshape(sum(toothF,2),n,3)-M.force_N;closure=max(abs(delta(:)));
assert(closure<1e-8,'R30F:Closure','Per-tooth force sum does not match total force.');
normToothN=sqrt(sum(toothN.^2,3));sumNormN=sum(normToothN,2);
S=r30f_spectrum([double(M.force_N),double(M.moment_Nmm(:,3))/1000,Fsingle,Fpair],dt);
fWheel=abs(M.parameters.installation.omegaW)/(2*pi);
fGear=abs(M.parameters.installation.omegaG)/(2*pi);
fMesh=fGear*M.parameters.p.gear.teeth;
if isempty(o.outputDir),o.outputDir=unique_folder(fullfile(root,['r30_academic_figures_' datestr(now,'yyyymmdd_HHMMSS')]));end
assert(~isfolder(o.outputDir),'R30F:Output','Output exists; choose a new outputDir.');mkdir(o.outputDir);
[ok,attr]=fileattrib(o.outputDir);assert(ok);out=attr.Name;
P=struct('kind','saved_force_thesis_figures_v2','sourceFile',sourceFile,'outputDir',out,...
 'forceDataChanged',false,'physicalCalibrationValidated',false,'singleTooth',o.singleTooth,...
 'pairTeeth',o.pairTeeth,'activeTeeth',ids(active),'candidateTeeth',ids,...
 'sampling',rmfield(S,{'frequency_Hz','amplitude'}),'nominalWheelRotation_Hz',fWheel,...
 'nominalGearRotation_Hz',fGear,'nominalMesh_Hz',fMesh,...
 'perToothForceClosure_N',closure,'pairOverlapDuration_s',nnz(norm1>0&norm2>0)*dt,...
 'pairPeak_N',max(normPair),'scalarSumPeak_N',max(upper),...
 'pairCancellationDefinition','eta=1-|F1+F2|/(|F1|+|F2|); NaN if denominator=0',...
 'pairMeaning','Vector superposition under prescribed rigid motion; no elastic load-redistribution model.',...
 'scope',M.scope,'figures',{{}},'plotTimeOrigin_s',double(M.window_s(1)),...
 'plotTimeMeaning','seconds since the saved common-window start; raw CSV retains absolute time',...
 'plotStyle','reference thesis: Chinese labels, Times numerals, inward ticks, blue/red thin lines, captions below panels');
write_csv(fullfile(out,'global_force_time.csv'),...
 {'time_s','Fx_N','Fy_N','Fz_N','F_resultant_N','Fn_vector_resultant_N','sum_per_tooth_normal_norm_N','Tc_Nm','nRemovingTeeth','nRemovingFaces'},...
 [t M.force_N vecnorm(M.force_N,2,2) vecnorm(M.normalVector_N,2,2) sumNormN M.moment_Nmm(:,3)/1000 double(M.nRemovingTeeth) double(M.nRemovingFaces)]);
write_csv(fullfile(out,'pair_vector_superposition.csv'),...
 {'time_s','F1x_N','F1y_N','F1z_N','F2x_N','F2y_N','F2z_N','Fpair_x_N','Fpair_y_N','Fpair_z_N','F1_norm_N','F2_norm_N','vector_sum_norm_N','sum_of_norms_N','cancellation_index'},...
 [t F1 F2 Fpair norm1 norm2 normPair upper eta]);
write_csv(fullfile(out,'fft_single_sided_amplitude.csv'),...
 {'frequency_Hz','global_Fx_N','global_Fy_N','global_Fz_N','global_Tc_Nm','single_Fx_N','single_Fy_N','single_Fz_N','pair_Fx_N','pair_Fy_N','pair_Fz_N'},...
 [S.frequency_Hz S.amplitude]);
headers={'time_s'};for k=1:numel(ids),headers{end+1}=sprintf('tooth_%02d_normal_resultant_N',ids(k));end
write_csv(fullfile(out,'per_tooth_normal_time.csv'),headers,[t normToothN]);
save(fullfile(out,'figure_numeric_data.mat'),'M','S','P','-v7');
write_text(fullfile(out,'FIGURE_NOTES_zh.txt'),notes(P));
if o.makeFigures
 P.figures=r30f_thesis_figures(M,S,o,out);
end
% Save report after native figure export is complete.
save(fullfile(out,'figure_report.mat'),'P','-v7');write_text(fullfile(out,'figure_report.json'),jsonencode(P));
fprintf('\n学术绘图完成：%s\n',out);
fprintf('N=%d, Fs=%.6f Hz, T=N*dt=%.6f s, FFT频率间隔=%.6f Hz。\n',n,S.sampleRate_Hz,S.recordDuration_s,S.binSpacing_Hz);
fprintf('单齿%d；双齿%d+%d，按矢量相加。未建立弹性载荷耦合模型。\n',o.singleTooth,o.pairTeeth);
fprintf('当前短窗只适合比较窗口内频谱；不能据此识别齿轮转频或固有频率。\n');
end

function out=unique_folder(base)
out=base;k=0;while isfolder(out),k=k+1;out=[base '_' num2str(k)];end
end
function write_csv(path,headers,data)
fid=fopen(path,'w');assert(fid>0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',strjoin(headers,','));
fmt=[repmat('%.17g,',1,size(data,2)-1) '%.17g\n'];fprintf(fid,fmt,data');
end
function write_text(path,value)
fid=fopen(path,'w','n','UTF-8');assert(fid>0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',value);
end
function text=notes(P)
text=sprintf(['图由保存的多齿M数据直接绘制；不平滑、不缩放，不改力模型。\n',...
 '采用论文风格：中文宋体、数字Times New Roman、黑框内向刻度、蓝红细线、图下分图说明；输出FIG、矢量PDF和400 dpi PNG。\n',...
 '图中时间从共同窗口起点计，CSV保留原始绝对时间。\n',...
 '全力为ΣF的模；纯法向合力为ΣN的模；Σ||N_i||是不同的标量上界。\n',...
 '双齿%d+%d逐分量矢量相加。相消指标η=1-||F1+F2||/(||F1||+||F2||)。\n',...
 '该图表示规定刚体运动下的载荷叠加，不表示弹性耦合或载荷重新分配。\n',...
 'FFT：原始均匀网格，先减普通均值，再乘对称Hann窗，单边幅值除以sum(w)，仅非DC/非Nyquist乘2。\n',...
 '无零填充，无平滑；幅值谱不是功率谱密度。DC不画出。\n',...
 'N=%d，Fs=%.9g Hz，T=N*dt=%.9g s，频率间隔Δf=%.9g Hz。Hann主瓣更宽，Δf不等于近邻频率可分辨能力。\n',...
 '名义砂轮转频%.9g Hz，名义啮合频率%.9g Hz；本工况二者重合，不能据峰位置区分来源。\n',...
 '当前窗口很短，不能解析约%.9g Hz的齿轮转频，也不能把谱峰直接认定为结构固有频率。\n',...
 '单齿载荷有进退场包络，短窗FFT不能替代完整稳定周期的阶次或频率分析。\n',...
 '未算历史的候选齿保存载荷为零，不表示已验证其剩余材料状态。\n',...
 '力系数未标定；名义轮径210 mm与当前求交轮径213.5764 mm仍待确认。\n'],...
 P.pairTeeth,P.sampling.sampleCount,P.sampling.sampleRate_Hz,P.sampling.recordDuration_s,P.sampling.binSpacing_Hz,...
 P.nominalWheelRotation_Hz,P.nominalMesh_Hz,P.nominalGearRotation_Hz);
end
