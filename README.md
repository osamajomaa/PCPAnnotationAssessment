This project aims to discover the hidden relationships between proteins in large databases and how they interact with each other. The database we use is Uniprot-GOA, which is basically a Uniprot but its protein sequences are annotation to the Gene Ontology. 

We construct networks entity-entity (EE nets) networks and entity-citation-entity networks (ECE nets). For netowkrs we have two types of entityes: proteins (PP and PCP) and GO terms (GG and GCG). In the EE nets, we have an edge between two entities (two proteins or tow GO terms) if they appear in the same paper. On the other hand, in the ECE nets, we have an edge between two entities if they appear in a paper and its citation. 

We visualize the four types of networks using matplotlib. We also visualize some statistic that help us analyze the networks and predict the proteins behavior. 
