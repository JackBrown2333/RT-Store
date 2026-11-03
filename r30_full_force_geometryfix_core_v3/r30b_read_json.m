function value=r30b_read_json(file)
% Explicit UTF-8, including when MATLAB runs under a Windows system locale.
fid=fopen(file,'rb');assert(fid>0,'R30B:Input','Cannot open JSON: %s',file);
cleanup=onCleanup(@()fclose(fid));bytes=fread(fid,Inf,'*uint8')';
value=jsondecode(native2unicode(bytes,'UTF-8'));
end
