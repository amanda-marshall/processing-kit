function structA = pk_mergeStructs(structA, structB)
% pk_mergeStructs: Merge two structures

% Returns a merged structure. If Structure B contains fields that also exist in
% structure A, fields from B will override A.

  f = fieldnames(structB);
  for i = 1:length(f)
    fieldname = f{i};
    % If both structs contain the same field and both are structs, deep merge them...
    if isfield(structA, fieldname) && isstruct(structA.(fieldname)) && isfield(structB, fieldname) && isstruct(structB.(fieldname))
      structA.(fieldname) = pk_mergeStructs(structA.(fieldname), structB.(fieldname));
    % ...otherwise, just overrride it
    else
      structA.(fieldname) = structB.(fieldname);
    end
  end

end
