function q_1=norm_q(q)
%Normalization
q_1=q;
sumq=sum(q);
for i=1:length(q)
    q_1(i)=q(i)/sumq;
end