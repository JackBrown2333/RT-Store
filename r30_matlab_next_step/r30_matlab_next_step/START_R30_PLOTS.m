% 打开本文件并按F5。自动加载本包，无需预先设置MATLAB当前目录。
r30LauncherDir=fileparts(mfilename('fullpath'));
assert(~isempty(r30LauncherDir),'R30:Launcher','请打开START_R30_PLOTS.m并按F5，不要只复制其内容。');
addpath(r30LauncherDir,'-begin');rehash;
r30Setup=r30f_setup_paths('plots');
fprintf('绘图入口：%s\n',which('run_r30_academic_figures'));
P=run_r30_academic_figures(r30Setup.resultDir);
