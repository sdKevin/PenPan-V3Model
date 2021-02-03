function [ day ] = Num_Days_Mon( x )
% 计算一年中每个月的天数 Num_Days_Mon(201501)
month = mod(x,100);
year = (x-month)./100;
if mod(month,2) && month<8
    day = 31;
elseif mod(month,2)==0&&month>=8
    day=31;
elseif mod(month,2)==0&&month<7&&month>3
    day=30;
elseif mod(month,2)&&month>8
    day=30;
else
    if mod(year,4)==0
        day=29;
    else day=28;
    end
end
end