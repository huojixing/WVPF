clear all
%This is the original matlab code of WVPF
%%
%DATA
[FileName,PathName] = uigetfile('*.xls*;*.lvm','OPEN EXCEL FILE');
Fxls=xlsread(strcat(PathName,FileName));
sF=size(Fxls);
N=100;%The number of particle filter
M_min=4;%l_min
M_max=271;%l_max
beta=0.3;
time=Fxls(:,1);%hour-min-sec
To=quling(Fxls(:,2));%Ourdoor temperature. Raw data should be processed
Ti=quling(Fxls(:,3));%Indoor temperature
l=length(time);
if l<M_max
    sheding=(max(Ti)-min(Ti))*beta;
else
    sheding=(max(Ti(1:M_max))-min(Ti(1:M_max)))*beta;
end
P=adP(Fxls(:,8));%Electrical power of air conditioner
Ph_additional=400;%Heat source
ch=-1;%Working mode
state_P=zeros(length(P),1);%On or off of air conditioner
for i=1:l
    if abs(P(i))>500
        P(i)=P(i)*ch;
        state_P(i)=1;
    else
        P(i)=0;
        state_P(i)=0;
    end
end
ts=rem(time,100);%Sec
tm=rem((time-ts)/100,100);%Min
th=rem(((time-ts)/100-tm)/100,100);%Hour
sTi=0.1;%Standard deviation of Ti
slgC=0.1/3;slgK=0.1/3;slga=0.01/3;%standard deviation of parameters
Ti_trace=Ti;%Predicted value
C=zeros(l,1);K=zeros(l,1);alpha=zeros(l,1);
%%
%Initialization
C(1)=2.88*10^6;K(1)=40;alpha(1)=4;
C_part=zeros(N,1);K_part=zeros(N,1);alpha_part=zeros(N,1);%Before Prediction particles
C_part_1=zeros(N,1);K_part_1=zeros(N,1);alpha_part_1=zeros(N,1);%Prediction particles
q=zeros(N,1);%weight
C_particle=10.^(log10(abs(C(1)))+slgC*2*(rand(N,1)-0.5));%Particle after resampling
K_particle=10.^(log10(abs(K(1)))+slgK*2*(rand(N,1)-0.5));
alpha_particle=10.^(log10(abs(alpha(1)))+slga*2*(rand(N,1)-0.5));
%%
%Proposed method WVPF
for k=1:N
    C_part(k)=C(1);
    K_part(k)=K(1);
    alpha_part(k)=alpha(1);
end
var_C_K_alpha=[slgC^2*ones(1,length(Ti));slgK^2*ones(1,length(Ti));slgK^2*ones(1,length(Ti))];%Variances in iteration
M_save=[];%Record the prediction window length
t_consumption=[];%Record the time consumption
for i=2:l-M_min+1
    tic
    if i>M_max
        sheding=(max(Ti(i-M_max:i))-min(Ti(i-M_max:i)))*beta;
    end
    %Compute prediction window length
    for M=M_min:min(M_max,l-i+1)
        if max(Ti(i:i+M-1))-min(Ti(i:i+M-1))>sheding
            break
        end
    end
    M_save=[M_save;M];
    h=ts(i)-ts(i-1);%Sampling period
    if h<0
        h=h+60;
    end
    %Prediction
    for k=1:N
        C_part_1(k)=10^(log10(C_part(k))+slgC*randn);
        K_part_1(k)=10^(log10(K_part(k))+slgK*randn);
        alpha_part_1(k)=10^(log10(alpha_part(k))+slga*randn);
    end
    Ti_pred=Ti_trace(i-1)*ones(N,1);q=ones(N,1);
    %Update
    for j=1:M
        for k=1:N
            %Fourth RK method
            K1=h*(-K_part_1(k)/C_part_1(k)*Ti_pred(k)+K_part_1(k)/C_part_1(k)*(To(i+j-2))+1/C_part_1(k)*alpha_part_1(k)*P(i+j-2)+1/C_part_1(k)*Ph_additional);
            K2=h*(-(K_part_1(k))/(C_part_1(k))*(Ti_pred(k)+0.5*K1)+(K_part_1(k))/(C_part_1(k))*((To(i+j-2)+To(i+j-1))/2)+1/(C_part_1(k))*alpha_part_1(k)*(P(i+j-2)+P(i+j-1))/2+1/C_part_1(k)*Ph_additional);
            K3=h*(-(K_part_1(k))/(C_part_1(k))*(Ti_pred(k)+0.5*K2)+(K_part_1(k))/(C_part_1(k))*((To(i+j-2)+To(i+j-1))/2)+1/(C_part_1(k))*alpha_part_1(k)*(P(i+j-2)+P(i+j-1))/2+1/C_part_1(k)*Ph_additional);
            K4=h*(-K_part_1(k)/C_part_1(k)*(Ti_pred(k)+K3)+K_part_1(k)/C_part_1(k)*(To(i+j-1))+1/C_part_1(k)*alpha_part_1(k)*P(i+j-1)+1/C_part_1(k)*Ph_additional);
            Ti_pred(k)=Ti_pred(k)+(K1+2*K2+2*K3+K4)/6;%predition of the true value of indoor temperature
            dTi=Ti(i+j-1)-Ti_pred(k);
            q(k)=q(k)*exp(-dTi^2/2/sTi^2);%Update the weight
        end
        if sum(q)==0%Avoid the normalization of weight is NULL
            q=ones(N,1);
            disp([i,j,M]);
        end
        q=norm_q(q);%Normalization in each predition instant
    end
    % Resample.
    for k=1:N 
        C_particle(k)=f_resample(C_part_1,q); 
        K_particle(k)=f_resample(K_part_1,q); 
        alpha_particle(k)=f_resample(alpha_part_1,q); 
    end
    % The particle filter estimate is the mean of the particles.
    C(i)=10^(mean(log10(abs(C_particle))));
    K(i)=10^(mean(log10(abs(K_particle))));
    alpha(i)=10^(mean(log10(abs(alpha_particle))));
    var_C_K_alpha(:,i)=[sum((log10(C_particle)-log10(C(i))).^2);sum((log10(K_particle)-log10(K(i))).^2);sum((log10(alpha_particle)-log10(alpha(i))).^2)]/N;
    %The predicted value is computed via 4th RK
    K1=h*(-K(i)/C(i)*Ti_trace(i-1)+K(i)/C(i)*(To(i-1))+1/C(i)*alpha(i)*P(i-1)+1/C(i)*Ph_additional);
    K2=h*(-K(i)/C(i)*(Ti_trace(i-1)+0.5*K1)+K(i)/C(i)*((To(i-1)+To(i))/2)+1/C(i)*alpha(i)*(P(i-1)+P(i))/2+1/C(i)*Ph_additional);
    K3=h*(-K(i)/C(i)*(Ti_trace(i-1)+0.5*K2)+K(i)/C(i)*((To(i-1)+To(i))/2)+1/C(i)*alpha(i)*(P(i-1)+P(i))/2+1/C(i)*Ph_additional);
    K4=h*(-K(i)/C(i)*(Ti_trace(i-1)+K3)+K(i)/C(i)*(To(i))+1/C(i)*alpha(i)*P(i)+1/C(i)*Ph_additional);
    Ti_trace(i)=Ti_trace(i-1)+(K1+2*K2+2*K3+K4)/6;
    %For iteration
    C_part=C_particle;K_part=K_particle;alpha_part=alpha_particle;
    t_c=toc;
    t_consumption=[t_consumption,t_c];
end
%²¹³ä
for i=l-M_min+2:l
    C(i)=C(l-M_min+1);K(i)=K(l-M_min+1);alpha(i)=alpha(l-M_min+1);
end
%%
%Predicted value of indoor temperature
Ti_com=Ti;
for i=2:l
    h=ts(i)-ts(i-1);
    if h<0
        h=h+60;
    end
    K1=h*(-K(i)/C(i)*Ti_com(i-1)+K(i)/C(i)*(To(i-1))+1/C(i)*alpha(i)*P(i-1)+1/C(i)*Ph_additional);
    K2=h*(-K(i)/C(i)*(Ti_com(i-1)+0.5*K1)+K(i)/C(i)*((To(i-1)+To(i))/2)+1/C(i)*alpha(i)*(P(i-1)+P(i))/2+1/C(i)*Ph_additional);
    K3=h*(-K(i)/C(i)*(Ti_com(i-1)+0.5*K2)+K(i)/C(i)*((To(i-1)+To(i))/2)+1/C(i)*alpha(i)*(P(i-1)+P(i))/2+1/C(i)*Ph_additional);
    K4=h*(-K(i)/C(i)*(Ti_com(i-1)+K3)+K(i)/C(i)*(To(i))+1/C(i)*alpha(i)*P(i)+1/C(i)*Ph_additional);
    Ti_com(i)=Ti_com(i-1)+(K1+2*K2+2*K3+K4)/6;
end

%%
%Output
l=3*720;
xx=0:5/3600:5/3600*(l-1);
figure(1)
plotyy(1:length(time),[Ti,To],1:length(time),ch*P)
figure(2)
plot(xx,C(1:l),'b')
xlim([0 l/720]);
xlabel('Time(h)');
ylabel('C(J/\circC)');
figure(3)
plot(xx,K(1:l),'b')
xlim([0 l/720]);
xlabel('Time(h)');
ylabel('K(W/\circC)');
figure(4)
plot(xx,alpha(1:l).*state_P(1:l))
figure(7)
plot([Ti_com,Ti])
RMSE=sqrt(sum((Ti-Ti_com).^2)/length(Ti))
