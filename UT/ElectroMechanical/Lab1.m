clear 
clc

groupA=[1701.725,1694.675,1705.955,1705.955,1703.135,1701.725,1701.725,1704.545,1701.725,1705.955,1701.725,1698.905,1697.495,1701.725,1703.135,1703.135,1697.495,1698.905,1704.545,1700.315,1701.725,1705.955,1704.545,1701.725,1701.725];
groupC=[1705.12,1704.00,1704.00,1704.00,1702.88,1705.12,1688.35,1710.71,1705.12,1705.12,1705.12,1705.12,1707.91,1702.32,1702.32,1705.12,1704.56,1705.12,1704.00,1705.68,1705.12,1705.12,1703.26,1706.98,1706.98];
groupD=[1702.362,1701.969,1714.961,1720.472,1713.780,1712.992,1713.780,1711.417,1716.142,1711.811,1714.961,1723.622,1721.654,1708.661,1712.205,1718.110,1718.110,1708.661,1708.268,1713.389,1722.441,1710.236,1716.535,1710.236,1712.992];
groupE=[1697.007,1702.559,1702.598,1702.598,1697.402,1697.441,1697.323,1691.811,1697.362,1708.465,1697.323,1691.889,1682.008,1687.559,1693.110,1682.008,1687.559,1687.559,1687.559,1687.559,1687.559,1693.110,1687.559,1687.559,1687.559];
groupF=[1702.346,1702,1704.76,1712.875,1712.875,1721.75,1712.875,1695.125,1695.125,1704,1704,1712.875,1686.25,1695.125,1704,1695.125,1704,1695.125,1704,1695.125,1704,1695.125,1695.125,1695.125,1704];

alldata=[groupA,groupC,groupD,groupE,groupF];

aavg=sum(groupA(1:25))/length(groupA);
cavg=sum(groupC(1:25))/length(groupC);
davg=sum(groupD(1:25))/length(groupD);
eavg=sum(groupE(1:25))/length(groupE);
favg=sum(groupF(1:25))/length(groupF);

a=0;c=0;d=0;e=0;f=0;

for i=1:25
    anew=((groupA(i)-aavg)^2)/24+a;
    a=anew;
end
stdA=sqrt(anew);


for i=1:25
    cnew=((groupC(i)-cavg)^2)/24+c;
    c=cnew;
end
stdC=sqrt(cnew);

for i=1:25
    dnew=((groupD(i)-davg)^2)/24+d;
    d=dnew;
end
stdD=sqrt(dnew);

for i=1:25
    enew=((groupE(i)-eavg)^2)/24+e;
    e=enew;
end
stdE=sqrt(enew);

for i=1:25
    fnew=((groupF(i)-favg)^2)/24+f;
    f=fnew;
end
stdF=sqrt(fnew);

xax=[1:25];
hold all
plot(xax,groupA,'-x','LineWidth',2,'markersize',15)
plot(xax,groupC,'LineWidth',2)
plot(xax,groupD,'LineWidth',2)
plot(xax,groupE,'LineWidth',2)
plot(xax,groupF,'LineWidth',2)
 

% figure

% histfit(groupA)
% figure
% histfit(groupC)
% figure
% histfit(groupD)
% figure
% histfit(groupE)
% figure
% histfit(groupF)
% figure
% histfit(alldata)


% 
% 
% figure
% std_dev = stdA;
% mean = aavg;
% x = mean-3*std_dev:.1:mean+3*std_dev;
% y = (std_dev*sqrt(2*pi))^-1*exp(-.5*((x-mean)/std_dev).^2);
% plot(x,y);
% 
% 






