% Open THIS new V3 entry and press F5 after copying the complete package.
r30V3Root=fileparts(mfilename('fullpath'));
assert(~isempty(r30V3Root),'R30FULL:V3Entry','请打开START_R30_FULL_FORCE_V3.m并按F5。');
assert(isfolder(fullfile(r30V3Root,'geometry_fix_v3'))&&isfile(fullfile(r30V3Root,'r30full_prepare_corefix.m')),...
 'R30FULL:V3Install','请完整解压V3包，保留子目录，不要只复制启动文件。');
addpath(r30V3Root,'-begin');clear r30full_prepare_corefix r30b_strips r30b_geometry;rehash;
run(fullfile(r30V3Root,'START_R30_FULL_FORCE.m'));
