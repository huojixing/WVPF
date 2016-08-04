function f=f_resample(s,q)%Multinomial Resampling
N=length(s);
u=rand; % uniform random number between 0 and 1    
qtempsum=0; 
for j=1:N 
    qtempsum=qtempsum+q(j); 
    if qtempsum>=u 
        f=s(j); 
        break; 
    end
end