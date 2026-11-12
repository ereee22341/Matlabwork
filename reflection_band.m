function b = reflection_band(f,R,f0,threshold)
% 提取包含 f0 的连续低反射频带，并插值求两侧阈值交点。
[~,i]=min(abs(f-f0)); assert(R(i)<threshold,'设计频率处的反射率未低于阈值。');
lo=i; hi=i;
while lo>1 && R(lo-1)<threshold, lo=lo-1; end
while hi<numel(f) && R(hi+1)<threshold, hi=hi+1; end
assert(lo>1 && hi<numel(f),'当前扫频范围未覆盖低反射频带的两个边界。');
fl=f(lo-1)+(threshold-R(lo-1))*(f(lo)-f(lo-1))/(R(lo)-R(lo-1));
fh=f(hi)+(threshold-R(hi))*(f(hi+1)-f(hi))/(R(hi+1)-R(hi));
b=[fl,fh,(fh-fl)/f0];
end
