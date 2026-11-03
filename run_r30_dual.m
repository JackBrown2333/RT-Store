function R30B=run_r30_dual(options)
% Unzip, cd to this folder, then: R30B=run_r30_dual;
% This recomputes independent selected-flank histories from t=0 in MATLAB.
if nargin<1,options=struct();end
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'matlab_recompute'));
if ~isfield(options,'outputDir'),options.outputDir=fullfile(root,'r30_dual_recomputed');end
R30B=r30b_run(fullfile(root,'parameters.json'),options);
end
