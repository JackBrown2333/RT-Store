function rows=r30b_force(rows,f)
% Kn is N/mm; reference lengths are mm. No granular model or legacy K1 units.
if isempty(rows),return;end
fields={'Vdot_mm3ps','Qprime_mm2ps','heq_mm','Fn_N','Ft_N','Fx_N','Fy_N','Fz_N',...
 'Nvec_x_N','Nvec_y_N','Nvec_z_N','Tvec_x_N','Tvec_y_N','Tvec_z_N','Mx_Nmm','My_Nmm','Mz_Nmm'};
for k=1:numel(fields),[rows.(fields{k})]=deal(0);end
for k=1:numel(rows)
 r=rows(k);assert(all(isfinite([r.b_mm r.lg_mm r.vc_mmps r.dt_s r.dV_mm3]))&&...
  min([r.b_mm r.lg_mm r.vc_mmps r.dt_s])>0&&r.dV_mm3>=0,'R30B:Strip','Invalid force inputs');
 r.Vdot_mm3ps=r.dV_mm3/r.dt_s;r.Qprime_mm2ps=r.Vdot_mm3ps/r.b_mm;r.heq_mm=r.Qprime_mm2ps/r.vc_mmps;
 r.Fn_N=f.Kn*r.b_mm*(r.heq_mm/f.h_ref)^f.m*(r.lg_mm/f.l_ref)^f.p+f.background_pressure*r.contact_area_mm2;
 r.Ft_N=f.mu*r.Fn_N;N=r.Fn_N*r.normal_response;T=r.Ft_N*r.tangent_response;F=N+T;
 M=r.Fn_N*r.normal_moment_response+r.Ft_N*r.tangent_moment_response;
 r.Fx_N=F(1);r.Fy_N=F(2);r.Fz_N=F(3);r.Nvec_x_N=N(1);r.Nvec_y_N=N(2);r.Nvec_z_N=N(3);
 r.Tvec_x_N=T(1);r.Tvec_y_N=T(2);r.Tvec_z_N=T(3);r.Mx_Nmm=M(1);r.My_Nmm=M(2);r.Mz_Nmm=M(3);rows(k)=r;
end
end
