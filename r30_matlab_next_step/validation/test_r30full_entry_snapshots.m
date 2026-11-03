function report=test_r30full_entry_snapshots(sourceMat,fixedCore)
% Real parameters/grid, independent nominal-stock snapshots near reported time.
% These are geometry/force regressions, NOT the actual continuous history.
addpath(fixedCore,'-begin');D=load(sourceMat,'M');M=D.M;
co=r30b_options(struct('ns',M.options.ns,'nz',M.options.nz,'order',M.options.order,...
 'dt_s',M.forceDt_s,'substeps',M.options.substeps,'forceStripWidth_mm',.05,'coefficients',M.coefficients,'plot',false));
times=[2.25 2.255 2.2573375 2.26 2.265 2.27 2.28 2.29];nonempty=0;maxClosure=0;maxForce=0;
for side=1:2
 base=r30b_init(M.parameters,co,[1 side]);before=r30b_heights(base);
 for t=times
  Q=base;ta=t-co.dt_s;
  for sub=1:co.substeps,[Q,~]=r30b_step(Q,ta+sub*(t-ta)/co.substeps,co);end
  after=r30b_heights(Q);v=sum(Q.weightsFine(:).*(before(:)-after(:)));
  rows=r30b_force(r30b_strips(Q,before,after,before,t,co.dt_s,co,false),co.coefficients);
  assigned=0;
  if ~isempty(rows)
   nonempty=nonempty+1;assigned=sum([rows.dV_mm3]);assert(all(isfinite([rows.Fn_N])));maxForce=max(maxForce,max([rows.Fn_N]));
  end
  maxClosure=max(maxClosure,abs(assigned-v));assert(abs(assigned-v)<1e-9);
 end
end
assert(nonempty>0,'R30FULL:EntryTest','所选入口姿态未发生去除，不能作为非零回归。');
report=struct('grid',[co.ns co.nz co.order],'tooth',1,'sides',[1 2],'times_s',times,...
 'independentNominalStockSnapshots',true,'actualContinuousHistory',false,'nonzeroForceCases',nonempty,...
 'maximumPartitionError_mm3',maxClosure,'maximumStripFn_N',maxForce,'actualUserFailureSnapshotAvailable',false);disp(report);
end
