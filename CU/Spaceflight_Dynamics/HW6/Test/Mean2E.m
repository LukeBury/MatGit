function [E]=Mean2E(M,e)

Ei=M-.2;
updateE=Ei-M;
E=Ei;
%compute E
while abs(updateE) > 1e-15
    updateE=(M-(E-e*sin(E)))/(1-e*cos(E));
    E =E + updateE;
end