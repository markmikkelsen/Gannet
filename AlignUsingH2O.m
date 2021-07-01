function [output MRS_struct] = AlignUsingH2O(input,MRS_struct)
%Align to water maximum (top-right plot)
A=size(input)
[number index] = max(abs(input),[],2);
index=index-A(1)/2;
for ii=1:A(2)
 output(:,ii)=circshift(input(:,ii),[-A 0]);
end

 MRS_struct.out.reject(:,MRS_struct.ii) = zeros(A(2),1);
end