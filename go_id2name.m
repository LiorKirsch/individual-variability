function go_cat_names = go_id2name(go_ids, GO)

go_cat_names = get(GO(go_ids).terms, 'name');

end