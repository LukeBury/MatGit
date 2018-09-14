clear
clc
clf

R=[4385.44085, 6578.15717, 8702.08397]
sep=[17.3333,18.2,19.6333]
reat=[15.8667,17.2667,18.1]
hold all
grid on
plot(R,sep, ':bv','LineWidth',2,'MarkerEdgeColor','k','markersize',15)
plot(R,reat, '--ko', 'LineWidth',2, 'MarkerEdgeColor','k','markersize',15)
set(gca,'xTick',[4000 5000 6000 7000 8000 9000])
set(gca,'YTick',[15.5 16 17 18 19 20])
