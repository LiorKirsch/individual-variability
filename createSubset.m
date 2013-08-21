load('humanOntologyObject.mat');

regionNames25 = {'Frontal Lobe', 'Parietal Lobe', 'Temporal Lobe' ,'Occipital Lobe' , 'Insula', 'Cingulate gyrus',...
'parahippocampal gyrus' , 'hippocampal formation' , 'Basal Ganglia','Basal Forebrain' ,'Claustrum', 'Amygdala'...
,'Thalamus','Subthalamus','Epithalamus','Hypothalamus','Mesencephalon','Cerebellum'...
 ,'Pons','Myelencephalon'}';

[a,b] = ismember(regionNames25, humanOntology.structureLabels(:,4) );

[subsetNodesIndex,c] = ismember( humanOntology.structureLabels(:,4), regionNames25);

save('subset2.mat', 'subsetNodesIndex');