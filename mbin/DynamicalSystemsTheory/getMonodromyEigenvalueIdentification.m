function [EV_pair1_1, EV_pair1_2, EV_pair2_1, EV_pair2_2, EV_equals1_1, EV_equals1_2] = getMonodromyEigenvalueIdentification(eigenvalues_row)
%%% Description
% Takes a row vector of the 6 eigenvalues of a monodromy matrix and
% identifies the 3 pairs. The pairs are [EV1, EV2], [EV3, EV4], and 
% [EV_equals1_1, EV_equals1_2]
%       
% ------------------------------------------------------------------------
%%% Inputs
% eigenvalues_row - [1x6] Row vector of monodromy matrix eigenvalues 
% ------------------------------------------------------------------------
%%% Outputs
% EV1          - [scalar] Eigenvalue 1
% EV2          - [scalar] Eigenvalue 2
% EV3          - [scalar] Eigenvalue 3
% EV4          - [scalar] Eigenvalue 4
% EV_equals1_1 - [scalar] First eigenvalue equal to 1
% EV_equals1_2 - [scalar] Second eigenvalue equal to 1
% ------------------------------------------------------------------------
% Created: 12/08/20
% Author : Luke Bury, luke.bury@colorado.edu
% ========================================================================
% % -------------------------------------------------
% %%% Setup
% % -------------------------------------------------
% %%% Rename for understanding of the reduction-style algorithm
remainingEigenvalues = eigenvalues_row;
% 
% % -------------------------------------------------
% %%% Find and remove the two cases of lambda = 1 (note, it's never equal
% %%% to 1 exactly
% % -------------------------------------------------
% %%% Complex eigenvalues
% logInd_nonComplex = ~imag(remainingEigenvalues);
% complexEigenvalues = remainingEigenvalues(~logInd_nonComplex);
% 
% %%% If there are only two eigenvalues with complex parts, they must be the
% %%% complex conjugate pair
% if length(complexEigenvalues)==2
%     EV_pair1_1 = complexEigenvalues(1);
%     EV_pair1_2 = complexEigenvalues(2);
%     
% %%% Else, if there are 4 eigenvalues with complex parts, then the 
% %%% periodicity pair must be buried in here
% elseif length(complexEigenvalues)==4
%     %%% First, let's look for the periodicity pair by differencing the real
%     %%% parts of these eigenvalues with 1
%     complexEigsMinus1 = abs(real(complexEigenvalues) - 1);
%     %%% Post-differencing, if the minimum value appears twice, let's call
%     %%% that our periodicity pair and we'll mark the other two as our
%     %%% complex conjugate pair
%     if sum(complexEigsMinus1==min(complexEigsMinus1))==2
%         complexEigsWithRealClosestTo1 = complexEigenvalues(complexEigsMinus1==min(complexEigsMinus1));
%         EV_equals1_1 = complexEigsWithRealClosestTo1(1);
%         EV_equals1_2 = complexEigsWithRealClosestTo1(2);
%         
%         complexEigPair = complexEigenvalues(complexEigsMinus1~=min(complexEigsMinus1));
%         EV_pair1_1 = complexEigPair(1);
%         EV_pair1_2 = complexEigPair(2);
%     %%% Pos-differencing, if the minimum value only appears once, let's
%     %%% call that one of our periodicity eigenvalues, then remove it from
%     %%% the stack, and find the next post-differencing minimum, which we'll
%     %%% call the other periodicitiy eigenvalue
%     elseif sum(complexEigsMinus1==min(complexEigsMinus1))==1
%         EV_equals1_1 = complexEigenvalues(complexEigsMinus1==min(complexEigsMinus1));
%         remainingComplexEigenvalues = complexEigenvalues(complexEigenvalues~=EV_equals1_1);
%         logIndex_otherMinima = abs(real(remainingComplexEigenvalues) - 1) == min(abs(real(remainingComplexEigenvalues) - 1));
%         if sum(logIndex_otherMinima)~=1
%             warning('Eigenvalue parsing problem: code black')
%         end
%         EV_equals1_2 = remainingComplexEigenvalues(logIndex_otherMinima);
%         logIndex_complexPair = ~logIndex_otherMinima;
%         complexPair = remainingComplexEigenvalues(logIndex_complexPair);
%         EV_pair1_1 = complexPair(1);
%         EV_pair1_2 = complexPair(2);
%     else
%         warning('Weird eigenvalue parsing error')
%     end
% %%% If there are no complex eigenvalues, then both pairs are inverse pairs    
% elseif isempty(complexEigenvalues)
%     d
% else
%     warning('Extra complex eigenvalue')
%     return
% end
% 
% %%% Should have either 2 or 4 remaining eigenvalues
% remainingEigenvalues = remainingEigenvalues(logInd_nonComplex);
% 
% %%% Pick out the inverse pair
% EV_pair2_1 = max(remainingEigenvalues);
% if length(EV_pair2_1) ~=1
%     warning('Too many max eigenvalues')
%     return
% end
% 
% %%% Should only have 1 or 3 remaining eigenvalues
% remainingEigenvalues = remainingEigenvalues(remainingEigenvalues~=EV_pair2_1);
% 
% %%% Complete finding the inverse pair which satisy EV3 = 1/EV4
% tempVec = abs(remainingEigenvalues - 1/EV_pair2_1);
% EV_pair2_2 = remainingEigenvalues(tempVec==min(tempVec));
% 
% %%% If the periodicity pair does not have a complex component, then it
% %%% hasn't been parsed yet
% if sum(logInd_nonComplex)>2
%     %%% The only two remaining eigenvalues should be the periodicity pair which
%     %%% equal with (with some numerical error)
%     remainingEigenvalues = remainingEigenvalues(remainingEigenvalues~=EV_pair2_2);
%     if length(remainingEigenvalues)~=2
%         warning('Eigenvalue parsing error')
%         return
%     end
%     EV_equals1_1 = remainingEigenvalues(1);
%     EV_equals1_2 = remainingEigenvalues(2);
% end


for jj = 1:2
    %%% Difference of all EVs with 1
    tempVec = abs(real(remainingEigenvalues) - 1);
    
    %%% Identify the two eigenvalues equal to 1
    % If they're equal so both get flagged by min()
    if sum((tempVec == min(tempVec))) == 2
        temp = remainingEigenvalues((tempVec == min(tempVec)));
        EV_equals1_1 = temp(1);
        EV_equals1_2 = temp(2);
        
    elseif jj == 1
        EV_equals1_1 = remainingEigenvalues((tempVec == min(tempVec)));
    elseif jj == 2
        EV_equals1_2 = remainingEigenvalues((tempVec == min(tempVec)));
    end
    
    %%% Only keep those that aren't the minimum
    remainingEigenvalues = remainingEigenvalues(~(tempVec == min(tempVec)));
    
    %%% Set this in case the first two minima are equal and both
    %%% get removed
    if isequal(size(remainingEigenvalues),[1,4])
        break
    end
end
% % % % % % 
% % % % % % % -------------------------------------------------
% % % % % % %%% Find EV1-EV4 for current eigenvalue set
% % % % % % % -------------------------------------------------
% % % % % % %%% Find lambda1 (we'll define it as the EV with the largest real
% % % % % % %%% component)
% % % % % % maxRealEV = remainingEigenvalues(real(remainingEigenvalues) == max(real(remainingEigenvalues)));
% % % % % % 
% % % % % % %%% Check if there were two maximums - it's possible that the largest 
% % % % % % %%% real component belonged to a complex conjugate pair who have the 
% % % % % % %%% same real component
% % % % % % if isequal(size(maxRealEV),[1,1]) %%% Single Value - EV_pair1_2 must follow real(EV_pair1_1) = 1/real(EV_pair1_2)
% % % % % %     %%% Set EV_pair1_1
% % % % % %     EV_pair1_1 = maxRealEV;
% % % % % % 
% % % % % %     %%% Remove EV_pair1_1 from remainingEigenvalues
% % % % % %     remainingEigenvalues = remainingEigenvalues(~(remainingEigenvalues == EV_pair1_1));
% % % % % % 
% % % % % %     %%% Finding EV_pair1_2 from real(EV_pair1_1) = 1/real(EV_pair1_2)
% % % % % % % %     tempVec = abs(real(EV_pair1_1) - 1./real(remainingEigenvalues));
% % % % % %     tempVec = abs(real(remainingEigenvalues) - 1/real(EV_pair1_1));
% % % % % %     EV_pair1_2 = remainingEigenvalues(tempVec == min(tempVec));
% % % % % % 
% % % % % %     %%% Remove EV_pair1_2 from remainingEigenvalues if it's [1x1]
% % % % % %     remainingEigenvalues = remainingEigenvalues(~(remainingEigenvalues == EV_pair1_2));
% % % % % % 
% % % % % % elseif isequal(size(maxRealEV),[1,2]) %%% EV_pair1_1 has complex conjugate 
% % % % % %     %%% Pick out the complex conjugate pair
% % % % % %     EV_pair1_1 = maxRealEV(1);
% % % % % %     EV_pair1_2 = maxRealEV(2);
% % % % % % 
% % % % % %     %%% Remove EV_pair1_2 from remainingEigenvalues
% % % % % %     remainingEigenvalues = remainingEigenvalues(~(real(remainingEigenvalues) == max(real(remainingEigenvalues))));
% % % % % % end
% % % % % % 
% % % % % % %%% Only two EVs left - they must be a pair
% % % % % % EV_pair2_1 = remainingEigenvalues(1);
% % % % % % EV_pair2_2 = remainingEigenvalues(2);


% -------------------------------------------------
%%% Find EV1-EV4 for current eigenvalue set
% -------------------------------------------------
%%% Find lambda1 (we'll define it as the EV with the largest real
%%% component)
maxRealEV = remainingEigenvalues(abs(real(remainingEigenvalues)) == max(abs(real(remainingEigenvalues))));

%%% Check if there were two maximums - it's possible that the largest 
%%% real component belonged to a complex conjugate pair who have the 
%%% same real component
if isequal(size(maxRealEV),[1,1]) %%% Single Value - EV_pair1_2 must follow real(EV_pair1_1) = 1/real(EV_pair1_2)
    %%% Set EV_pair1_1
    EV_pair1_1 = maxRealEV;

    %%% Remove EV_pair1_1 from remainingEigenvalues
    remainingEigenvalues = remainingEigenvalues(~(remainingEigenvalues == EV_pair1_1));

    %%% Finding EV_pair1_2 from real(EV_pair1_1) = 1/real(EV_pair1_2)
% %     tempVec = abs(real(EV_pair1_1) - 1./real(remainingEigenvalues));
    tempVec = abs(real(remainingEigenvalues) - 1/real(EV_pair1_1));
    EV_pair1_2 = remainingEigenvalues(tempVec == min(tempVec));

    %%% Remove EV_pair1_2 from remainingEigenvalues if it's [1x1]
    remainingEigenvalues = remainingEigenvalues(~(remainingEigenvalues == EV_pair1_2));

elseif isequal(size(maxRealEV),[1,2]) %%% EV_pair1_1 has complex conjugate 
    %%% Pick out the complex conjugate pair
    EV_pair1_1 = maxRealEV(1);
    EV_pair1_2 = maxRealEV(2);

    %%% Remove EV_pair1_2 from remainingEigenvalues
    remainingEigenvalues = remainingEigenvalues(~(abs(real(remainingEigenvalues)) == max(abs(real(remainingEigenvalues)))));
end

%%% Only two EVs left - they must be a pair
EV_pair2_1 = remainingEigenvalues(1);
EV_pair2_2 = remainingEigenvalues(2);

end % function