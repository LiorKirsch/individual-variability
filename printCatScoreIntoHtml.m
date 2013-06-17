function printCatScoreIntoHtml(regionLabel, scores, catIds, catLabels, numOfgenesInCatAndRegion, title, htmlFileName )

    brainRelated = load('/home/lab/noalis/work3/for_uri_feb_2012/new/old_data/images_go_genes_mat_brain_screened', 'go_cat_names', 'cat_ids');
    
    if ~exist('htmlFileName','var')
        htmlFileName = 'output.html';
    end
    
    if ~exist('catLabels','var')
        GO = geneont('file', 'gene_ontology.obo.txt');
        catLabels = go_id2name(catIds, GO); 
    end

    fid = fopen(fullfile('www',htmlFileName), 'w');
    fprintf('writing output to file %s\n', htmlFileName);
    fprintf(fid,'<html><head><title>%s</title></head><body><table>\n', title);
    
    fprintf(fid,'<tr>');
    for i = 1:length(regionLabel)
        fprintf(fid,'<th> %s </th>', regionLabel{i});
    end
    fprintf(fid,'</tr>\n');

    maxIter = min(100, size(scores,1));
    for i = 1:maxIter
        fprintf(fid,'<tr>');
        for j = 1:length(regionLabel)
            link = sprintf('http://amigo.geneontology.org/cgi-bin/amigo/term_details?term=GO:%07d',catIds(i,j));
            
            style = '';
            if ismember(catIds(i,j), brainRelated.cat_ids)
                style = [style ,' background-color: #cfcfcf;'];
            end
            
            if scores(i,j) > 0.05
                style = [style ,'visibility:hidden;'];
            end
            
            fprintf(fid,'<td style="%s"> <small>%d</small> <a href=%s>%s</a> (%d) </td>', style, scores(i,j),link, catLabels{j}{i}, numOfgenesInCatAndRegion(i,j) );
        end
        fprintf(fid,'</tr>\n');
        printPercentCounter(j, size(scores,1));
    end
    
    %fprintf(fid, '</table></body></html>');
    fclose(fid);
end