function my_table = make_table(planetNames, u_of_planet, values)
% make a table of planet properties 
%
% INPUTS
%   planetNames   : cell array of planet names (rows)
%   u_of_planet   : cell array, u_of_planet{pp} = wind speeds for planet pp
%   values        : numeric array or cell array values(pp,i) = final value at wind i for planet pp (e.g., wave heights)
%
% OUTPUT
%   my_table      : 
%                   rows    = planets
%                   columns = wind speeds
%                   entries = values (NaN if not available) (e.g., wave heights)


    all_winds = unique([u_of_planet{:}]);
    all_winds = sort(all_winds);

    varNames = strcat('U_', strrep(string(all_winds),'.','p')); % Note: 2.2 = 2p2


    my_table = array2table(NaN(numel(planetNames), numel(all_winds)), 'RowNames', planetNames(:), 'VariableNames', varNames );
    for pp = 1:numel(planetNames)
        winds = u_of_planet{pp};

        for ii = 1:numel(winds)
            colName = strcat('U_', strrep(string(winds(ii)),'.','p'));
            my_table{planetNames{pp}, colName} = values(pp,ii);
        end
    end
end
