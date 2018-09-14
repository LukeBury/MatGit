function [dvs] = dJC_2_dv(dJCs)
%%% Inputs:
% 1) Jacobi constant difference (eg, JC2-JC1) [nx1]
%%% Outputs:
% 1) delta-v required to achieve jacobia constant difference
% ========================================================================
dvs = zeros(size(dJCs));
for kk = 1:length(dJCs)
    if dJCs(kk) < 0
        dvs(kk) = 0;
        warning('No dV necessary to theoretically reach position')
    else
        dvs(kk) = sqrt(dJCs(kk));
    end
end

end