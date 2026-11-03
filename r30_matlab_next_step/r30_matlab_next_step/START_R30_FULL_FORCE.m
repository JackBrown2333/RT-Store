% 解压本包全部内容到已有r30_matlab_next_step目录，打开本文件按F5。
% 全30齿、60独立齿面，从名义毛坯开始，保留完整轴向进给及真实齿号相位。
r30FullRoot=fileparts(mfilename('fullpath'));assert(~isempty(r30FullRoot));addpath(r30FullRoot,'-begin');rehash;
assert(isfile(fullfile(r30FullRoot,'r30f_setup_paths.m')),'R30FULL:Install','请把本包全部内容放进已有r30_matlab_next_step目录。');
clear run_r30_full_force r30full_collect r30full_read r30full_figures test_r30full_stream r30full_prepare_corefix test_r30full_centroid;
R30_FULL_FORCE_RESULT=struct('completed',false,'fullGearComputed',false);
r30FullSetup=r30f_setup_paths('direction');
r30FullPrep=fullfile(r30FullSetup.resultDir,'full_process_preparation_v1');
r30FullPlanFile=fullfile(r30FullPrep,'preflight','full_trial_plan.mat');
assert(isfile(r30FullPlanFile)&&isfile(fullfile(r30FullPrep,'resolution','resolution_metrics.json')),...
 'R30FULL:Preparation','缺少本机完整行程预检/收敛报告，请先运行START_R30_FULL_PREP。');
r30FullCheck=jsondecode(fileread(fullfile(r30FullPrep,'resolution','resolution_metrics.json')));
assert(strcmp(r30FullCheck.status,'COMPLETED_DIFFERENCES_REPORTED')&&r30FullCheck.maximumBaselineForceDifference_N<1e-7&&r30FullCheck.maximumBaselineMomentDifference_Nmm<1e-5);
r30FullOriginalCore=r30FullSetup.coreDir;
r30FullSetup.coreDir=r30full_prepare_corefix(r30FullOriginalCore,fullfile(r30FullSetup.projectDir,'r30_full_force_geometryfix_core_v3'));
r30FullOut=fullfile(r30FullSetup.projectDir,'r30_full_force_30T_160x320_b005_geometryfix_v3');
if ~isfolder(r30FullOut),mkdir(r30FullOut);end
r30FullZip=fullfile(r30FullOut,'r30_full_force_summary.zip');if isfile(r30FullZip),delete(r30FullZip);end
r30FullSource=fullfile(r30FullSetup.resultDir,'multitooth_pilot.mat');
% Cache software self-test only when every tested compute source is identical.
addpath(fullfile(r30FullRoot,'validation'),'-begin');
r30FullSignature='';
for r30FullName={'run_r30_full_force.m','r30full_collect.m','r30full_read.m','validation/test_r30full_stream.m','validation/test_r30full_centroid.m'}
 r30FullSignature=[r30FullSignature fileread(fullfile(r30FullRoot,r30FullName{1}))];
end
r30FullCoreFiles=dir(fullfile(r30FullSetup.coreDir,'*.m'));
for r30FullJ=1:numel(r30FullCoreFiles)
 r30FullSignature=[r30FullSignature r30FullCoreFiles(r30FullJ).name fileread(fullfile(r30FullSetup.coreDir,r30FullCoreFiles(r30FullJ).name))];
end
r30FullCache=fullfile(r30FullOut,'software_selftest_cache.mat');r30FullTestNeeded=true;
if isfile(r30FullCache)
 r30FullCacheData=load(r30FullCache,'r30FullSavedSignature','r30FullTestReport');
 r30FullTestNeeded=~isequal(r30FullCacheData.r30FullSavedSignature,r30FullSignature);
end
if r30FullTestNeeded
 r30FullTestDir=tempname(r30FullOut);mkdir(r30FullTestDir);
 fprintf('先执行小规模软件自检（合成工况），通过后自动转入30齿真实工况的全程试算。\n');
 r30FullCentroidReport=test_r30full_centroid(r30FullOriginalCore,r30FullSetup.coreDir,r30FullSource);
 r30FullTestLog=evalc('r30FullTestReport=test_r30full_stream(r30FullSource,r30FullSetup.coreDir,r30FullTestDir);');
 r30FullLogFile=fopen(fullfile(r30FullTestDir,'selftest_log.txt'),'w','n','UTF-8');
 if r30FullLogFile>0,fprintf(r30FullLogFile,'%s',r30FullTestLog);fclose(r30FullLogFile);end
 r30FullSavedSignature=r30FullSignature;save(r30FullCache,'r30FullSavedSignature','r30FullTestReport','-v7');
 fprintf('软件自检通过：分块载荷、最终材料状态及续算一致。\n');
end
try
 r30FullRun=run_r30_full_force(r30FullSource,r30FullSetup.coreDir,r30FullOut,r30FullPlanFile);
 assert(r30FullRun.completed&&r30FullRun.fullGearComputed,'R30FULL:Incomplete','全30齿尚未完成。');
 r30FullFiles={'full_run_report.json','full_run_result.mat','full_force_overview.csv','final_face_statistics.csv',...
  'peak_window.mat','central_segment_fft.mat','central_segment_fft_0_1000Hz.csv'};
 for r30FullFigure={'full_force_overview','peak_multitooth_forces','contact_and_remaining_stock','central_segment_fft'}
  r30FullFiles{end+1}=fullfile('figures',[r30FullFigure{1} '.png']);
 end
 r30FullData=load(fullfile(r30FullOut,'full_run_result.mat'),'stats');
 [~,r30FullWorstRemain]=max(r30FullData.stats(:,5));[~,r30FullWorstOvercut]=max(r30FullData.stats(:,6));
 r30FullShow=unique([17 r30FullData.stats(r30FullWorstRemain,1) r30FullData.stats(r30FullWorstOvercut,1)]);
 for r30FullTooth=r30FullShow,r30FullFiles{end+1}=fullfile('figures',sprintf('final_tooth_%02d.png',r30FullTooth));end
 r30FullExists=cellfun(@(name)isfile(fullfile(r30FullOut,name)),r30FullFiles);r30FullFiles=r30FullFiles(r30FullExists);
 zip(r30FullZip,r30FullFiles,r30FullOut);assert(isfile(r30FullZip));
 R30_FULL_FORCE_RESULT=r30FullRun;R30_FULL_FORCE_RESULT.resultZipFile=r30FullZip;
catch r30FullError
 fprintf(2,'全程任务尚未完成。已提交的材料检查点保留，再运行同一入口可继续。\n');
 fprintf(2,'实际错误：%s [%s]\n',r30FullError.message,r30FullError.identifier);
 if any(strcmp(r30FullError.identifier,{'R30B:NormalLayer','R30:NormalLayer','R30:RootBracket'}))
  fprintf(2,'本次为法向材料求交错误；需检查几何，不能截平目标面。\n');
 end
 if isfile(fullfile(r30FullOut,'last_force_error.mat'))
  r30FullErrorFiles={'last_force_error.mat','full_run_plan.mat'};
  r30FullErrorZip=fullfile(r30FullOut,'r30_full_force_error.zip');zip(r30FullErrorZip,r30FullErrorFiles,r30FullOut);
  fprintf(2,'请回传这个错误诊断包：\n%s\n',r30FullErrorZip);
 end
 rethrow(r30FullError);
end
fprintf('\n全程试算完成，请上传汇总ZIP：\n%s\n',r30FullZip);
fprintf('原始载荷块、全部最终材料状态、全部FIG/PDF保留在本机。\n');
if ispc,try,winopen(r30FullOut);catch,fprintf('结果目录：%s\n',r30FullOut);end;end
