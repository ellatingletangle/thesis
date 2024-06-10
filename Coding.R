rm(list = ls())
setwd("~/Documents/Master P&T/Thesis/Bundestag") 

library("rJava")
library("rDNA")
library("statnet")
library("igraph")
library("cluster")
library("remotes")
library("network")
library("ggplot2")
library("dplyr")
library("reshape2")

set.seed(420)
dna_init(jarfile="dna-2.0-beta25.jar")
#Open the DNA graphical user interface:
#dna_gui("Codes_new.dna")

#Connect to database
data <- dna_connection("Codes_new.dna", verbose = T)

###### Overview #######
concept_frequency <- dna_getAttributes(data, statementType= "DNA Statement", variable = "concept")
organization_frequency <- dna_getAttributes(data, statementType= "DNA Statement", variable = "organization")
references_frequency <- dna_getAttributes(data, statementType= "Scientific Reference ", variable = "Evidence")

sum(concept_frequency$frequency[concept_frequency$type == "Secondary Aspects"])
sum(concept_frequency$frequency[concept_frequency$type == "Policy core beliefs "])

#Plot frequencies 
ggplot(concept_frequency, aes(frequency, reorder(value, +frequency))) +
  geom_col(width=0.5, mapping=aes(fill = type))+
  scale_fill_manual(name="Type", values = c("magenta","orange"), 
                  labels=c("Policy Core Beliefs", "Secondary Aspects")) + 
  labs(x= "Frequency", y= "Concept") + 
  theme_bw()+
  theme(legend.position = c(0.7, 0.2))

ggplot(organization_frequency, aes(frequency, reorder(value,+frequency))) +
  geom_col(width=0.5, mapping=aes(fill = type)) +
  scale_fill_manual(name="Type",  
                    values = c("green", "red", "magenta", "blue","yellow"), 
                    labels = c("Civil Society Organization", "Economic Actor", "Labor Union", "Political Actor", "Research Institution")) + 
  labs(x= "Frequency", y= "Organization") + 
  theme_bw()+
  theme(legend.position = c(0.7, 0.3))

###### Network Analysis #####
#General Agreement over all concepts 
dna_barplot(data, of = "concept", fontSize = 10, colors=T, barWidth = 0.5, truncate = 100,  excludeValues = list(concept = ""))

#Only agreement of organizations with the Concept Cannabis Legalization
dna_barplot(data,
            of = "organization",
            fontSize = 5,
            excludeValues = list(concept = "Proposed policies lack supporting evidence "),
            invertValues = TRUE,
            colors = T
            )

dna_barplot(data,
            of = "concept",
            fontSize = 5,
            excludeValues = list(organization = "SPD"),
            invertValues = TRUE,
            colors = T
)
#Exclude: All concepts with zero disagreement.  
excluded_concept <- list(concept = c( "Consumption is increasing ", "Harm reduction approach is best method",  "Adolescent cannabis use is harmful and should be prevented.", "Drug laws major source for societal problems ", "Impurities in cannabis products are dangerous ", "High THC-concentration in cannabis is dangerous ", "Law enforcement is overwhelmed ", "Cannabis major source for societal problems ", "Increased addiction potential of cannabis is a problem ", "Societal circumstances cause drug use"))

### Coalition network ###
congruence <- dna_network(data,
                          networkType = "onemode",
                          statementType = "DNA Statement",
                          variable1 = "organization",
                          variable2 = "concept",
                          qualifier = "agreement",
                          qualifierAggregation = "congruence",
                          duplicates = "document", 
                          excludeValues = excluded_concept) 

dna_plotNetwork(congruence, layout = "fr", edge_weight = T,  edge_size_range= c(0.5,5),node_attribute = "type" ,node_colors = "manual", custom_colors = c("green", "red", "magenta", "blue","yellow"),
                node_label = T, theme = "graph",  threshold = 5, show_legend=T) 


######### Cluster Analysis ########
congruence_graph <- dna_toIgraph(congruence, weighted = TRUE)
# Perform community detection using the fast greedy algorithm
cluster <- cluster_fast_greedy(congruence_graph)

E(congruence_graph)$weight
V(congruence_graph)$color<- organization_frequency$color
V(congruence_graph)$type <- organization_frequency$type
# Define edge line type based on weight; if weight is less than 5, set to 0 (invisible), otherwise set to 1 (solid)
edge_lty <- ifelse(E(congruence_graph)$weight < 5, 0, 1)

set.seed(420)
plot(cluster, congruence_graph, 
     col= V(congruence_graph)$color, 
     vertex.frame.color =  V(congruence_graph)$color, 
     vertex.label.family = "sans",
     vertex.label.color= "black",
     vertex.label.cex = 0.7,
     mark.col=c("pink", "lightgreen"), 
     mark.border= c( "pink", "lightgreen"), 
     edge.width = E(congruence_graph)$weight/10, 
     edge.lty =  edge_lty,
     layout= layout.fruchterman.reingold(congruence_graph)
     )
legend(0.1,  
       legend= c("Economic Actor","Opposition Party", "Civil Society Organization","Labor Union", "Research Institution", "Governing Party"), 
       pch = 21,
       col = unique(V(congruence_graph)$color),
       pt.bg = unique(V(congruence_graph)$color),
       pt.cex = 2,
       cex = 0.7,
       x.intersp = 0.5,
       y.intersp = 0.5,
       bty = "n"
       )

########## Belief system of communities ###########
affil <- dna_network(data,
                     networkType = "twomode",
                     statementType = "DNA Statement",
                     variable1 = "organization",
                     variable2 = "concept",
                     qualifier = "agreement",
                     qualifierAggregation = "combine",
                     duplicates = "document", 
                     excludeValues = excluded_concept) 
                   
affiliations <- network(affil, bipartite = T)
# Extract membership information from a clustering object
memberships<-as.vector(membership(cluster))
# Create a data frame with node names and their corresponding community memberships
nodes_df <- data.frame(node = V(congruence_graph)$name, community = memberships)
# Convert the affiliation matrix to a data frame
affiliation_df <- as.data.frame(affil)
affiliation_df$actor <- rownames(affil)
# Merge affiliation data frame with node community data frame
affiliation_community_df <- merge(affiliation_df, nodes_df, by.x = "actor", by.y = "node")
# Calculate the mean agreement for each community, handling NA values appropriately
community_affiliation <- affiliation_community_df %>%
  group_by(community) %>%
  mutate(across(where(is.numeric), ~replace(., . == 0, NA))) %>% 
  summarise(across(`Agreement with current policies`:`Support current or additional crime penalties for cannabis use. `, ~ mean(.x, na.rm = TRUE)))
# Reshape the data frame to a long format suitable for plotting  
df_long <- melt(community_affiliation, id.vars = 'community', variable.name = 'Variable'
                , value.name = 'Value')
df_long$community<- as.factor(df_long$community)
df_long$Value <- df_long$Value - 1.5 
df_long$Value <- df_long$Value * - 2
# Create a bar plot
ggplot(df_long, aes(Variable, Value)) +
  geom_bar(stat= "identity", position = "stack", mapping = aes(fill = community)) +
  scale_fill_manual(name="Coalition",  values = c("#FC7878","lightgreen"), 
                    labels = c( "Prohibiton coalition","Legalisation coalition"
                    )) + 
  scale_y_continuous(breaks = seq(-1, 1, by=0.5), labels = c("Disagree","Partially disagree", "Neutral", "Partially  agree" ,"Agree")) + 
  coord_flip() +
  labs(y = "Agreement", x = "Concepts")+
  geom_hline(yintercept=0)+
  theme_bw()+
  theme(legend.position = "top")

###### Conflict Analysis ########
congruence_graph <- delete_edges(congruence_graph, E(congruence_graph)[weight < 5])

### Calculating within-group density for each cluster ###
intra_group_densities <- sapply(membership(cluster), function(x) {
  subg <- induced_subgraph(congruence_graph, which(membership(cluster) == x))
  possible_edges <- vcount(subg) * (vcount(subg) - 1) / 2  #55 for community 1,136 for 2
  actual_edges <- ecount(subg) #53 for community 1, 136 for 2
  density <- actual_edges / possible_edges
  return(density)
})

### Calculating between-group density for each pair of clusters ###
# Get edge IDs between the two groups
cluster1_nodes <- which(membership(cluster) == 1)
cluster2_nodes <- which(membership(cluster) == 2)

# Create all combinations of pairs between the two clusters
pairs <- expand.grid(cluster1_nodes, cluster2_nodes)
# Convert to a vector, with each pair followed by the next
vp <- as.vector(t(pairs))
cross_edges <- get.edge.ids(congruence_graph, vp)
# Count the number of edges between groups (ignoring weights)
num_cross_edges <- length(cross_edges[!cross_edges==0])
possible_cross_edges <- length(cross_edges)
# Calculate density as the ratio of actual to possible connections
between_group_density <- num_cross_edges / possible_cross_edges
print(paste("Between-group Density:", between_group_density))

### Density of conflict network ###
conflict <- dna_network(data,
                            networkType = "onemode",
                            statementType = "DNA Statement",
                            variable1 = "organization",
                            variable2 = "concept",
                            qualifier = "agreement",
                            qualifierAggregation = "conflict",
                            duplicates = "document", 
                            excludeValues = excluded_concept) 

dna_plotNetwork(conflict, layout = "fr", edge_weight = T,  edge_size_range= c(0.5,5), node_attribute = "type" ,node_colors = "manual", custom_colors = c("green", "red", "magenta", "blue","yellow"), node_label = T , threshol =2, theme = "graph", show_legend=T) 

conflict_network <- dna_toIgraph(conflict, weighted = TRUE)
E(conflict_network)$weight
ecount(conflict_network)

conflict_network <- delete_edges(conflict_network, E(conflict_network)[weight < 2])
# Calculate density
conflict_network_density <- ecount(conflict_network) / (28*27)/2


######## Scientific References ########
science <- dna_getStatements( data, statementType= "Scientific Reference ") 

docs <- dna_getDocuments(data)
actor <- orgaization_frequency$value
orga <- c("CDU/CSU",  "AfD", "SPD", "CDU/CSU", "CDU/CSU", "SPD",  "Die Linke", "Die Grünen", "SPD", "FDP", "Die Grünen"  ,"SPD", "FDP", "CDU/CSU", "SPD",  "CDU/CSU", "CDU/CSU", "AfD", "CDU/CSU", "SPD", "CDU/CSU", "SPD", "Die Grünen" ,  "Die Grünen" ,"Die Linke", "SPD" , "FDP", "AfD",  "DHV", "DPolG", "BPC" , "DHS"  , "Schildower Kreis", "CSCD",  "BAH", "NRV", "ZIS", "BvCW", "DGKJP", "DGS", "DR","GHN",  "BAK","BPtK","CDR" , "BVKJ" ,"Akzept","VCA" , "DGPPN", "ABDA") 
docs$organization <- orga
docs <- rename(docs, "documentId" = "id")

# Merge using left_join to keep all records from 'science'
evidence_df <- science %>%
  left_join(docs, by = "documentId")

#Count the number of times each actor made a reference
count_df <- evidence_df %>%
  group_by(organization, Evidence) %>%
  summarise(Count = n(), .groups = 'drop')  
print(count_df)

count_df$coalition <- c(1,1,2,1,1,1,2,1,1,1,2,2,2,2,2,2,2,2,2,2,2)

ggplot(count_df, aes(reorder(organization, coalition), Count)) +
  geom_col(width=0.5, mapping=aes(fill = Evidence), position = "dodge")+
  scale_fill_manual(name="Evidence Type", values = c("#BC7BCC","#AAB8DC"), 
                    labels=c( "Evidence opposing legalization", "Evidence supporting Legalization")) + 
  labs(x= "Organization", y= "Count") + 
  geom_vline(xintercept = 7.5) +
  theme_bw()+
  theme(legend.position = "top")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  

    