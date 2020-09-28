function structB = pk_mergeStructs(structA, structB)
% pk_mergeStructs: Merge two structures

% Returns a merged structure. If Structure B contains fields that also exist in
% structure A, fields from B will override A.

  f = fieldnames(structA);
  for i = 1:length(f)
     structB.(f{i}) = structA.(f{i})
  end

end;
