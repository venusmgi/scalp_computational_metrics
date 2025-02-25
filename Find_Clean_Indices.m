function [indices] = Find_Clean_Indices(binary_sequence, num_consecutive_points)


count = 0;
j = 1;
indices = [];
for i = 1: length(binary_sequence)
    if binary_sequence(i) ==1 
        count = count +1;
        if mod(count, num_consecutive_points)==0
            indices(end+1)= i - num_consecutive_points + 1;  %%figure out how youcan add 330,1073,1597,1822 it has something to do with diff of 4, but then, 1081 is calculated because you added it
            j = j+1;
            count = 0;
        end
    else 
        count = 0;
    end
end
end