# Fitness-position relationship (after TN-seq)

From experience and literature we know that genes near replication origin are present with an higher number of copy w.r.t. genes which are more distant; this indirectly affect the fitness value for those genes, pumping them up a little bit more than others. This behaviour can be graphically seen in the genomic_position-fitness plot, where the the distribution of points in regions near replication orgin (0) assume a 'smile-shape' behaviour.

The objective of this project is to correct this smiley behaviour since genes fitness should not be dependent on the genomic position. After the correction some observations are made.

In this project we consider data of genomic positions and fitness of several genes in _Streptococcus pneumoniae_.
Data were collected from this article: https://www.nature.com/articles/nmeth.1377

We have fitness values and genomic positions in 2 different files, respectively Tn_seq_fitness_data_Opijnen_et_al_2009.txt and GCF_000006885.1_ASM688v1_genomic_olt.txt.
