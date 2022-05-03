load('/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/GEM/personalized_GEM/GroupFBA/Supple_File4_hcGEM.mat')

hcGEM = addReaction(hcGEM, 'ATP_Hydrolysis','reactionFormula','m01371c + m02040c => m01285c + m02039c + m02751c','reversible',false,'subSystem','Transport reactions','checkDuplicate',true,'lowerBound',0,'upperBound',1000);

essentialTasks = parseTaskList('/home/anoop/Desktop/COVID_Omics/GEM_4_OriginalGrps/metabolicTasks_Essential.xlsx');
checkTasks(hcGEM, [], true, false, false, essentialTasks);


hcGEM = simplifyModel(hcGEM);
hcGEM.b = hcGEM.b(:,1);
hcGEM = setParam(hcGEM,'obj','ATP_Hydrolysis',1);

[NUM,STR] = xlsread('/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/GEM/personalized_GEM/Bounds/BoundsHC.xlsx');
cRxns = STR(2:end,1); 
cLB = NUM(:,1);
cUB = NUM(:,2);

hcGEM = setParam(hcGEM,'lb',cRxns,cLB);
hcGEM = setParam(hcGEM,'ub',cRxns,cUB);

sol_hcGEM = solveLP(hcGEM);
writematrix(sol_hcGEM.x,'/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/GEM/personalized_GEM/GroupFBA/Fba_hcGEM.txt')
writecell(hcGEM.rxns,'/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/GEM/personalized_GEM/GroupFBA/Reaction_hcGEM.txt')

clear

%%%%%%%%%%%
load('/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/GEM/personalized_GEM/GroupFBA/Suppl_File3_ecGEM.mat')

ecGEM = addReaction(ecGEM, 'ATP_Hydrolysis','reactionFormula','m01371c + m02040c => m01285c + m02039c + m02751c','reversible',false,'subSystem','Transport reactions','checkDuplicate',true,'lowerBound',0,'upperBound',1000);

essentialTasks = parseTaskList('/home/anoop/Desktop/COVID_Omics/GEM_4_OriginalGrps/metabolicTasks_Essential.xlsx');
checkTasks(ecGEM, [], true, false, false, essentialTasks);


ecGEM = simplifyModel(ecGEM);
ecGEM.b = ecGEM.b(:,1);
ecGEM = setParam(ecGEM,'obj','ATP_Hydrolysis',1);

[NUM,STR] = xlsread('/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/GEM/personalized_GEM/Bounds/BoundsEC.xlsx');
cRxns = STR(2:end,1); 
cLB = NUM(:,1);
cUB = NUM(:,2);

ecGEM = setParam(ecGEM,'lb',cRxns,cLB);
ecGEM = setParam(ecGEM,'ub',cRxns,cUB);

sol_ecGEM = solveLP(ecGEM);
writematrix(sol_ecGEM.x,'/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/GEM/personalized_GEM/GroupFBA/Fba_ecGEM.txt')
writecell(ecGEM.rxns,'/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/GEM/personalized_GEM/GroupFBA/Reaction_ecGEM.txt')

clear
%%%%%%
load('/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/GEM/personalized_GEM/GroupFBA/Suppl_File5_artGEM.mat')

artGEM = addReaction(artGEM, 'ATP_Hydrolysis','reactionFormula','m01371c + m02040c => m01285c + m02039c + m02751c','reversible',false,'subSystem','Transport reactions','checkDuplicate',true,'lowerBound',0,'upperBound',1000);

essentialTasks = parseTaskList('/home/anoop/Desktop/COVID_Omics/GEM_4_OriginalGrps/metabolicTasks_Essential.xlsx');
checkTasks(artGEM, [], true, false, false, essentialTasks);


artGEM = simplifyModel(artGEM);
artGEM.b = artGEM.b(:,1);
artGEM = setParam(artGEM,'obj','ATP_Hydrolysis',1);

[NUM,STR] = xlsread('/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/GEM/personalized_GEM/Bounds/BoundsART.xlsx');
cRxns = STR(2:end,1); 
cLB = NUM(:,1);
cUB = NUM(:,2);

artGEM = setParam(artGEM,'lb',cRxns,cLB);
artGEM = setParam(artGEM,'ub',cRxns,cUB);

sol_artGEM = solveLP(artGEM);
writematrix(sol_artGEM.x,'/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/GEM/personalized_GEM/GroupFBA/Fba_artGEM.txt')
writecell(artGEM.rxns,'/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/GEM/personalized_GEM/GroupFBA/Reaction_artGEM.txt')

clear
