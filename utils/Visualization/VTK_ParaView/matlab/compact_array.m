function new_array=compact_array(old_array)
% this function compacts an old array to a new
% array by sorting it and removing the repeated
% entries

temp = sort(old_array);
for i = 1: length(temp)
  if (i == 1)
    new_array = [temp(1)];
  else 
    if temp(i) ~= temp(i-1)
       new_array = [new_array,temp(i)];
    end
  end
end

