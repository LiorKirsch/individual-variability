function [sorted_C,sorted_ia,sorted_ib] =intersect_sort_by_B(A,B)

    [C,ia,ib] = intersect(A,B);
    
    [sorted_ib, sort_ind] = sort(ib);
    sorted_ia = ia(sort_ind);
    sorted_C = C(sort_ind);
    
end