function writeParFile(varargin)

if nargin > 1
    if rem(numel(varargin),2)
        error('Number of inputs must be an even number, array of named values ...');
    end
else
    error('Not enough input arguments ...');
end

for i = 1:2:numel(varargin)
    if ~ischar(varargin{i}) || ~isa(varargin{i+1},'double')
        error('Inputs must be an array of named values ...');
    end
end

File = fopen('data.par','w');
if ~File
    error('File can not be opened ...');
end

for i = 1:2:numel(varargin)
    
    WriteName(varargin{i});
    fwrite(File,varargin{i+1},'double');
end

fclose(File);

    function WriteName(name)
        for j = 1:length(name)
            fwrite(File,name(j),'char');
        end
        for j = length(name)+1:80
            fwrite(File,' ','*char');
        end
    end

end