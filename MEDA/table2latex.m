function table2latex(T, fileName, formatspec, begend)

% table2latex(T, 'table.tex', '%.2f', 0);
% Function to convert the table structure that is output from the MEDA toolbox to a latex markup for publication.
% Works with the Octave table structure, but assumes that the input is an Octave table structure if Octave is being used.
% Michael Sorochan Armstrong, 2022-12-19
%
% Create exception for edge case where Matlab is being used, but there is an Octave table structure.

isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;

if ~isOctave
    T2.source = table2cell(T(:,'Source'));
    T2.var = T.Properties.VariableNames;
    T2.mat = table2array(T(:,2:end));
    T = T2;
end

if nargin < 4
    begend = 0;
end

fid = fopen(fileName, 'w');

if begend == 0 || begend == 1
    fprintf(fid,'\\begin{tabular}{');
    fprintf(fid,repmat('l',[1 size(T.var,2) + 1]));
    fprintf(fid,'}\n');
end

if begend == 0 || begend == 1 || begend == 2 || begend == -1
    for ii = 1:size(T.mat,1)+1
        for jj = 1:size(T.mat,2)+1
            if jj == 1 && ii == 1
                fprintf(fid, ' &');
                fprintf(fid, ' ');
            elseif jj == 1
                fprintf(fid,strcat(T.source{ii-1,1},' &'));
                fprintf(fid,' ');
            elseif ii == 1 && jj == size(T.mat,2)+1
                fprintf(fid, T.var{1,jj});
                fprintf(fid,' ');
            elseif ii == 1
                fprintf(fid, strcat(T.var{1,jj},' &'));
                fprintf(fid,' ');
            elseif jj == size(T.mat,2)+1
                fprintf(fid, num2str(T.mat(ii-1,jj-1),formatspec));
                fprintf(fid,' ');
            else
                fprintf(fid,strcat(num2str(T.mat(ii-1,jj-1),formatspec),' &'));
                fprintf(fid,' ');
            end
        end

        if ii == 1
            fprintf(fid,'\\\\ \n');
            fprintf(fid,' \\hline \n');
        else
            fprintf(fid,'\\\\ \n');
        end
    end
end

if begend == 0 || begend == 2
    fprintf(fid,'\\end{tabular} \n');
end

fclose(fid);
