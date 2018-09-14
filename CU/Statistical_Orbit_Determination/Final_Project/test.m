clear
clc

syms x y z xs ys zs

r = [x y z];
rs = [xs ys zs];

range = sqrt(dot((r-rs),(r-rs)))
simplify(range)


jacobian(range,r)