function [Lia,Locb] = ismember_sort_by_B(A,B)

[Lia,Locb] = ismember(A,B);
Locb_relevent = sort( Locb(Lia) );

Locb = zeros(size(Locb));
Locb(Lia) = Locb_relevent;

end