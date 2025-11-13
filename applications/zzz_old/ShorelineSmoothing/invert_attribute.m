function STRUCT = invert_attribute(STRUCT)

fields = fieldnames(STRUCT);
for i = 1:numel(fields)
    if isnumeric(STRUCT.(fields{i})) && ismatrix(STRUCT.(fields{i}))
        STRUCT.(fields{i}) = transpose(STRUCT.(fields{i}));
    end
end

end
