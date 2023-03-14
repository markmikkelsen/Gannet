function [output, MRS_struct] = AlignUsingH2O(input, MRS_struct)
% Align to residual water magnitude

A = size(input);
[~, index] = max(abs(input),[],1);
index = index - A(1)/2;

for ii = 1:A(2)
    output(:,ii) = circshift(input(:,ii), -index(ii), 1); %#ok<AGROW> 
end

MRS_struct.out.reject{MRS_struct.ii} = zeros(1,A(2));

fprintf('\n');

end
