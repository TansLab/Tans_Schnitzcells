function yesno = ...
    inlist(mylist,element)

yesno = ~isempty(find(mylist==element));