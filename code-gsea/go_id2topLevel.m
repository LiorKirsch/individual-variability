function go_cat_topLevelHier = go_id2topLevel(go_ids, GO)

    go_cat_topLevelHier = get(GO(go_ids).terms, 'ontology');

end