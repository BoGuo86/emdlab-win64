
function y = getlist(commanname,index)

n = length(index);
y = cell(1,n);

for i = 1:n
    y{i} = [commanname,num2str(index(i))];
end

end