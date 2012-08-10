function result = isdiagonal (a)
[i,j] = find(a);
if ~isempty(i)
    result = all(i == j);
else
    result = true;
end
end