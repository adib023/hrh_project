function [varargout] =read(fname)
fid = fopen ( fname, "r");

i=1;

while ~feof(fid)
    tline = fgetl(fid);
    data_cell {i} = str2num (tline);
    i= i +1;
end

flag = 1;
new_cell_num = 1;

for n = 1 : size ( data_cell,2)-1
    
    if size(cell2mat(data_cell(n))) == [0 0]
        
        data = zeros (size(cell2mat(data_cell(n+1))));
        
        continue
        
    elseif size(cell2mat(data_cell(n))) ~= size(cell2mat(data_cell(n+1)))
        
        data(flag,:) = cell2mat(data_cell(n));
        new_cell {new_cell_num} = data;
        new_cell_num = new_cell_num + 1;
        flag = 1;
        
    else
        data(flag,:) = cell2mat(data_cell(n));
        flag = flag+1;
        
    end
    
end

for n =  1 : nargout
    varargout{n} = cell2mat(new_cell(n));
end

fclose ( fid );

end
