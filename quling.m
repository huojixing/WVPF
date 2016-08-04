function T=quling(T_1)
%Process the raw data of temperature
T=T_1;
z_location=[];
for i=1:length(T)
    if T(i)==0
        z_location=[z_location,i];
        if i==length(T)
            dT=T(z_location(1)-1)-T(z_location(1)-2);
            for j=1:length(z_location)
                T(z_location(j))=T(z_location(1)-1)+dT*j;
            end
            z_location=[];
        elseif T(i+1)~=0
            if z_location(1)==1
                dT=T(z_location(length(z_location))+2)-T(z_location(length(z_location))+1);
                for j=1:length(z_location)
                    T(z_location(z_location-j+1))=T(z_location(length(z_location))+1)-dT*j;
                end
            else
                dT=(T(z_location(length(z_location))+1)-T(z_location(1)-1))/(length(z_location)+1);
                for j=1:length(z_location)
                    T(z_location(j))=T(z_location(1)-1)+dT*j;
                end
            end
            z_location=[];
        end
    end
    if i>1&&abs(T(i)-T(i-1))>1
        T(i)=T(i-1);
    end
end
