function[array_unique] = SearchSame(array)
m = size(array,2);
count = 1;
array_unique = [];
for i = 1 : m
    for j = i+1 : m
         if(array(:,i) == array(:,j))
            array(:,j) = [];
            m = m-1;
            count = count + 1;
         end
         if(j >= m)
            break;
         end

    end
    if(i == m+1)
        break;
    end
   unique_temp = [array(:,i);count];
   array_unique = [array_unique,unique_temp];
   count = 1;
end  
end