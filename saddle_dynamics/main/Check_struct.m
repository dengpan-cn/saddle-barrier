function Idx = Check_struct(Set, simInfo)
% Check_struct.m - Check whether a configuration already exists in a set.
%
% Author: Yuchen Xie
% Copyright (c) 2024-2026 Yuchen Xie.
%
% Inputs:
% Set:          array of archived configurations
% simInfo:      configuration to be compared against Set
%
% Output:
% Idx:          index of the matched configuration in Set; 0 if not found

n = length(Set);
Idx = 0;
i = 1;
while i <= n
    flag = jampel_compareConf(Set(i), simInfo);
    if flag
        break
    else
        i = i + 1;
    end
end
if i <= n
    Idx = i;
end

end
