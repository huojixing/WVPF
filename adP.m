%Process the raw data of electrcal power
function P=adP(P)
for i=1:length(P)
    if P(i)==-1&&i==1
        P(i)=P(2);
    end
    if P(i)==-1&&i>1
        P(i)=P(i-1);
    end
end