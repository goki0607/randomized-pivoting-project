function test()
    %files = {'beaconfd'; 'cari'; 'farm'; 'itest6'; ...
    %    'klein2'; 'nsic1'; 'nsic2'; 'osa-07'; 'osa-14'; ...
    %    'osa-30'; 'rosen1'; 'rosen2'; 'rosen8'; 'rosen10'; ...
    %    'sc205'; 'scfxm1'; 'sctap2'; 'sctap3'};
    %files = {'boeing2'};
    %files = {'beaconfd'};
    files = dir('../tests/*.mps');
    in = files(1).folder;
    rows = {'danzig'; 'bland'; 'steepest'; 'randedge'; ...
        'randfacet'; 'clarkson'};
    cols = {'t'; 'f'; 'its'; 'flag'; 'fails'};
    csv = '.csv';
    out = '../results/';
    success = {};
    failure = {};
    for i=1:length(files)
        %file = char(fullfile(in,strcat(files(i),mps)));
        name = files(i).name;
        file = strcat(in,strcat('/',name));
        p = mpsread(file);
        arr = zeros(length(rows),length(cols));
        flag = 0;
        j = 1;
        while j <= length(rows)
            disp(j)
            try
                [~,t,f,its,flag] = evallp(p,j-1);
                arr(j,1:4) = [t f its flag];
                lastwarn('');
                [warnMsg, ~] = lastwarn;
                if ~isempty(warnMsg)
                    arr(5) = arr(5) + 1;
                    if arr(5) <= 10
                        continue
                    end
                end
                j = j + 1;
            catch
                flag = 1;
                break
            end
        end
        if flag == 1
            failure{end+1} = [name, rows(j)]; %#ok<AGROW>
            continue
        end
        success{end+1} = name; %#ok<AGROW>
        table = array2table(arr,'VariableNames',cols,'RowNames',rows);
        name = name(1:end-4);
        outfile = char(fullfile(out,strcat(name,csv)));
        disp(outfile)
        writetable(table, outfile, ...
            'Delimiter',',','QuoteStrings',true)
    end 
    table = cell2table(success);
    outfile = char(fullfile(out,strcat('success',csv)));
    writetable(table, outfile, ...
         'Delimiter',',','QuoteStrings',true)
    table = cell2table(failure);
    outfile = char(fullfile(out,strcat('failure',csv)));
    writetable(table, outfile, ...
         'Delimiter',',','QuoteStrings',true)
end