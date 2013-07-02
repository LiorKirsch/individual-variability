function printCatScoreIntoHtml(regionLabel, scores, scoresSeperate, catIds, catLabels, go_gene_mat, geneNames, numOfgenesInCatAndRegion, title, htmlFileName ,normalizedMeanGeneScore)

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
    collapseScript = '<script type="text/javascript" src="CollapsibleLists.js"></script>';
    fprintf(fid,'<html><head><title>%s</title>%s</head><body><table>\n', title, collapseScript);
    
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
            
            go_to_gene = logical(go_gene_mat{j});
            genesInCategory = go_to_gene(i,:);
            geneListHtml = createGeneString(geneNames, genesInCategory, normalizedMeanGeneScore(:,j));
            scoresIndHtml = createScoresSeperate( scoresSeperate(i,:,j) );
            fprintf(fid,'<td style="%s"> <ul class="collapsibleList">  <li> <small>%d</small> %s <a href=%s>%s</a> (%d) %s </li>  </ul></td>', style, scores(i,j),scoresIndHtml, link, catLabels{j}{i}, numOfgenesInCatAndRegion(i,j) ,geneListHtml);
        end
        fprintf(fid,'</tr>\n');
        printPercentCounter(j, size(scores,1));
    end
    
    fprintf(fid,'</table>\n<script type="text/javascript">     CollapsibleLists.apply();   </script>');
    fprintf(fid, '</body></html>');
    fclose(fid);
end

function scoresListHtml = createScoresSeperate(scoresSeperate)
            
            scoresListHtml = '';
            for m = 1:length(scoresSeperate)
                scoresListHtml = sprintf('%s <li>%g</li>', scoresListHtml, scoresSeperate(m) );
            end
            scoresListHtml = sprintf('<small> <ol>%s</ol> </small>' ,scoresListHtml);
end


function geneListHtml = createGeneString(geneNames, genesInCategory, normalizedMeanGeneScore)
            
            genesInCat = geneNames( genesInCategory );
            geneScores = normalizedMeanGeneScore( genesInCategory );
            [sortedScores, sortInd] = sort(geneScores);
            genesInCat = genesInCat(sortInd);
            geneListHtml = '';
            for m = 1:length(genesInCat)
                geneListHtml = sprintf('%s <li> %g <a href="http://www.genecards.org/cgi-bin/carddisp.pl?gene=%s">%s </a> </li> ', geneListHtml, sortedScores(m), genesInCat{m}, genesInCat{m});
            end
            geneListHtml = sprintf('<ul> %s </ul>' ,geneListHtml);
end