function s = mergestructs(s1,s2)
% s = mergestructs(s1,s2)
% 
% Merge two structures, s1 and s2, into a single structure, s

s = s1;
fields2 = fieldnames(s2);
for ii = 1:numel(fields2)
    s.(fields2{ii}) = s2.(fields2{ii});
end

end
