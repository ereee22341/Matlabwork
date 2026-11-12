function b = reflection_band(f,R,f0,threshold)
% Continuous band containing f0; interpolate both threshold crossings.
[~,i]=min(abs(f-f0)); assert(R(i)<threshold);
lo=i; hi=i;
while lo>1 && R(lo-1)<threshold, lo=lo-1; end
while hi<numel(f) && R(hi+1)<threshold, hi=hi+1; end
assert(lo>1 && hi<numel(f),'Band edges not bracketed');
fl=f(lo-1)+(threshold-R(lo-1))*(f(lo)-f(lo-1))/(R(lo)-R(lo-1));
fh=f(hi)+(threshold-R(hi))*(f(hi+1)-f(hi))/(R(hi+1)-R(hi));
b=[fl,fh,(fh-fl)/f0];
end
