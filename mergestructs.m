function s = mergestructs(s1,s2)
% s = mergestructs(s1,s2)
% 
% Merge two structures, s1 and s2, into a single structure, s

s = cell2struct([struct2cell(s1); struct2cell(s2)], [fieldnames(s1); fieldnames(s2)]);

end
