function r30b_outputs(R,folder)
write_curve(fullfile(folder,'combined_for_supplied_faces.csv'),R.poses);
for f=1:size(R.faces,1),write_curve(fullfile(folder,sprintf('tooth%d_side%d.csv',R.faces(f,:))),R.facePoses{f});end
meta=rmfield(R,{'poses','facePoses'});fid=fopen(fullfile(folder,'run_metadata.json'),'w','n','UTF-8');assert(fid>0);
fprintf(fid,'%s',jsonencode(meta));fclose(fid);
if ~R.options.plot,return;end
try
visible='off';if R.options.visible,visible='on';end
fig=figure('Visible',visible,'Color','w');t=[R.poses.t_s];
subplot(4,1,1);hold on;labels=cell(size(R.faces,1)+1,1);
for f=1:size(R.faces,1),plot(t,[R.facePoses{f}.Fn_pure_resultant_N],'LineWidth',1.1);labels{f}=sprintf('Tooth %d side %d',R.faces(f,:));end
plot(t,[R.poses.Fn_pure_resultant_N],'k-','LineWidth',1.5);labels{end}='Combined pure-normal vector norm';hold off;
legend(labels,'Location','best');grid on;ylabel('Pure normal load / N');title('Recomputed from t=0; candidate coefficients');
subplot(4,1,2);F=vertcat(R.poses.force_N);plot(t,F,'LineWidth',1.1);grid on;legend('Fx','Fy','Fz');ylabel('Force / N');
subplot(4,1,3);plot(t,[R.poses.F_reference_compression_N],t,[R.poses.Fn_magnitude_sum_N]);
grid on;legend('Full-force reference projection','Normal magnitude sum');ylabel('Load / N');
subplot(4,1,4);plot(t,[R.poses.Tc_Nm]);grid on;ylabel('C-axis moment / Nm');xlabel('Time / s');
saveas(fig,fullfile(folder,'dual_force_history.png'));savefig(fig,fullfile(folder,'dual_force_history.fig'));if ~R.options.visible,close(fig);end
catch e
 fid=fopen(fullfile(folder,'graphics_error.txt'),'w');
 if fid>0,fprintf(fid,'%s\n',e.message);fclose(fid);end
 warning('R30B:Graphics','Numerical MAT/CSV results were saved; plotting failed: %s',e.message);
end
end

function write_curve(file,S)
names={'t_start_s','t_s','dt_s','dV_mm3','Vdot_mm3ps','area_mm2','nForceCells',...
 'Fx_N','Fy_N','Fz_N','Nx_N','Ny_N','Nz_N','Fn_pure_resultant_N','Fn_magnitude_sum_N',...
 'F_reference_compression_N','N_reference_compression_N','Mx_Nmm','My_Nmm','Mz_Nmm','Tc_Nm',...
 'b_sum_mm','heq_widthWeighted_mm','lg_widthWeighted_mm','ae_pass_peak_mm','dh_step_peak_mm'};
fid=fopen(file,'w');assert(fid>0);cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s',names{1});for k=2:numel(names),fprintf(fid,',%s',names{k});end;fprintf(fid,'\n');
for k=1:numel(S)
 s=S(k);v=[s.t_start_s,s.t_s,s.dt_s,s.dV_mm3,s.Vdot_mm3ps,s.area_mm2,s.nForceCells,...
  s.force_N,s.normalVector_N,s.Fn_pure_resultant_N,s.Fn_magnitude_sum_N,...
  s.F_reference_compression_N,s.N_reference_compression_N,s.moment_Nmm,s.Tc_Nm,...
  s.b_sum_mm,s.heq_widthWeighted_mm,s.lg_widthWeighted_mm,s.ae_pass_peak_mm,s.dh_step_peak_mm];
 fprintf(fid,'%.16g',v(1));fprintf(fid,',%.16g',v(2:end));fprintf(fid,'\n');
end
end
