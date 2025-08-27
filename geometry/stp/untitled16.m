% Read STEP file
fid = fopen('result.step', 'r');
lines = textscan(fid, '%s', 'Delimiter', '\n', 'Whitespace', '');
fclose(fid);
lines = lines{1};

% Merge multi-line STEP entities
statements = {};
buffer = '';

for i = 1:length(lines)
    line = strtrim(lines{i});
    if isempty(line)
        continue
    end
    buffer = [buffer ' ' line];
    if endsWith(line, ';')
        statements{end+1} = strtrim(buffer); %#ok<AGROW>
        buffer = '';
    end
end

% Write new STEP file with one-line entities
fid = fopen('result1.step', 'w');
for i = 1:length(statements)
    fprintf(fid, '%s\n', statements{i});
end
fclose(fid);

disp('STEP file rewritten with single-line entities.');
