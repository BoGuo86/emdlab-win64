function str = emdlab_flib_getTextArray(x)

str = "[";
for i = 1:length(x)
    str = str + sprintf('%f',x(i)) + ",";
end
str = char(str);
str(end) = [];
str = string(str);
str = str + "]";

end