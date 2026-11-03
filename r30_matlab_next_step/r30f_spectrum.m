function S=r30f_spectrum(x,dt)
% Single-sided AMPLITUDE spectrum, not PSD. No toolbox or zero padding.
% Columns are independent real signals. Ordinary mean removed, then symmetric
% Hann window; coherent-gain correction is 1/sum(window). DC/Nyquist not doubled.
assert(isnumeric(x)&&isreal(x)&&all(isfinite(x(:))),'R30F:Signal','Signal must be finite and real.');
assert(isscalar(dt)&&isfinite(dt)&&dt>0,'R30F:Dt','dt must be positive.');
if isvector(x),x=x(:);end
n=size(x,1);assert(n>=4,'R30F:Length','At least four samples required.');
w=.5-.5*cos(2*pi*(0:n-1)'/(n-1));meanValue=mean(x,1);
Y=fft(bsxfun(@times,bsxfun(@minus,x,meanValue),w),[],1);
q=floor(n/2)+1;amplitude=2*abs(Y(1:q,:))/sum(w);amplitude(1,:)=amplitude(1,:)/2;
if mod(n,2)==0,amplitude(end,:)=amplitude(end,:)/2;end
S=struct('frequency_Hz',(0:q-1)'/(n*dt),'amplitude',amplitude,...
    'sampleRate_Hz',1/dt,'recordDuration_s',n*dt,'binSpacing_Hz',1/(n*dt),...
    'sampleCount',n,'fftLength',n,'window','symmetric Hann',...
    'coherentGain',sum(w)/n,'removedMean',meanValue,...
    'meaning','single-sided amplitude, not PSD; no zero padding; short finite record',...
    'note','Bin spacing is not the ability to resolve close spectral lines; Hann broadens peaks.');
end
