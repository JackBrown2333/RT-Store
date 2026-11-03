function checks=test_r30f_numerics(sourceFile,out)
% Numerical tests; no MATLAB graphics backend is required.
checks=struct();n=8000;dt=1/8000;t=(0:n-1)'*dt;
s=r30f_spectrum(3.2*sin(2*pi*100*t+.4)+10,dt);
[v,k]=max(s.amplitude(2:end));k=k+1;
assert(s.frequency_Hz(k)==100&&abs(v-3.2)<1e-6);checks.sinusoid_amplitude_error=abs(v-3.2);
s=r30f_spectrum(7*ones(n,1),dt);assert(max(s.amplitude)==0);checks.constant_removed=true;
s=r30f_spectrum(2*(-1).^(0:n-1)',dt);assert(abs(s.amplitude(end)-2)<1e-12);checks.nyquist_not_doubled=true;
n=7999;x=1.7*cos(2*pi*73*(0:n-1)'/n);s=r30f_spectrum(x,dt);
assert(numel(s.frequency_Hz)==4000&&abs(s.amplitude(74)-1.7)<1e-6);checks.odd_length=true;
s=r30f_spectrum([x 2*x],dt);assert(max(abs(s.amplitude(:,2)-2*s.amplitude(:,1)))<1e-12);checks.independent_columns=true;
P=run_r30_academic_figures(sourceFile,struct('makeFigures',false,'visible',false,'outputDir',out));
assert(P.sampling.sampleCount==8800&&P.sampling.fftLength==8800);
assert(abs(P.sampling.binSpacing_Hz-1/.11)<1e-9&&P.perToothForceClosure_N<1e-8);
checks.actual_M_8800_samples=true;checks.no_zero_padding=true;checks.force_sum_closure_N=P.perToothForceClosure_N;
checks.note='Numerics tested with GNU Octave; native MATLAB FIG/PDF export requires running the script in MATLAB.';
fid=fopen(fullfile(out,'numerical_test_report.json'),'w');fprintf(fid,'%s',jsonencode(checks));fclose(fid);
disp(checks);
end
