function res = bandpass_ifft( a, f1, f2 )
% 带通滤波，通带为[f1, f2]，f1<f2
% 借助傅里叶变换达到快速滤波
n=length(a);
if mod(n,2) == 0
    f = [ 1:floor(n/2), floor(n/2)-1:-1:1]/n;
else  %
    f = [ 1:floor(n/2), floor(n/2):-1:1]/n;
end
ind = f<f1|f>f2; % 获得阻带各个频率的下标
F = fft(a);
F2 = F(2:end);  
F2(ind)=0; % 把阻带响应设置为0
F(2:end)=F2; 
res = real( ifft(F) ); %逆变换回到时域
end