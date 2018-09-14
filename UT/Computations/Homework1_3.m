clear
clc

%Defining the variables
i = 2:2:10;
l_1 = 0;
l_2 = 0;

%Calculating L1 and L2
x=10.^-i;
l_1=sin(x)./x;
l_2=(1-cos(x))./x.^2;

%Output
for i=1:5
    a=[l_1' l_2'];
    fprintf('The limit of L1 is %.1g \t The limit of L2 is %.1g \n', a(i,:))
end