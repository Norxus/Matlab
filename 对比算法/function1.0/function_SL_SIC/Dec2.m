function [des] = Dec2(Dec,M,NT)
Dec = Dec - 1;
des = [];
while(Dec > 0)
    temp = mod(Dec,M);
    des=[temp;des];
    Dec=(Dec-temp)/M;
end
if(length(des)<NT)
    for i = 1 : (NT-length(des))
        des = [0;des];
    end
end
end