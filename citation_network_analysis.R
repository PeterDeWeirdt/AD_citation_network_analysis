#PETER DEWEIRDT
#Last Update: 9/27 - Generalizing Code
#CITATION NETWORK ANALYSIS

#Install packages we need
library(Matrix)
library(RSQLite)
library(igraph)
#packages for wordcloud analysis
library("tm")
library("SnowballC")
library("wordcloud")
library(plyr)
library(randomcoloR)
library(xtable)
library(doSNOW)
library(tcltk)
library(pbapply)
library(Hmisc)
library(NMI)

#Global Parameters
field = "title, abstract" #where do we want to draw words from?
db_path = "path to SQL database"
#Connect to our database 
db = dbConnect(SQLite(), dbname=db_path)
db = dbConnect(db)
#Get a read for what's in our database
dbListTables(db)
dbListFields(db,"paper_citation") 

#Network Creation
network = create_network()
ecount(network) 
vcount(network) 
paper_ids = clusters_table$V1
network = set_vertex_attr(network, "paper_ids", value = V(network))
### If we just want Alzheimer's Papers - uncomment and replace network with alz_network###
#query = "SELECT id FROM papers WHERE degree = 0"
#alz_paper_ids = dbGetQuery(db, query)
#alz_vertices = which(network$paper_ids %in% alz_paper_ids$id)
#alz_network = induced.subgraph(whole_network, alz_vertices, impl = "auto")
#alz_network = induced_subgraph(alz_network, which(components(alz_network)$membership == 1))
#alz_network = simplify(alz_network)
### Otherwise ###
network = induced_subgraph(network, which(components(network)$membership == 1))
network = simplify(network)
attrs = get.vertex.attribute(network)
paper_ids = attrs$paper_ids
years = get_paper_years(paper_ids)
initial_probabilities = get_initial_probs(as.numeric(years))
cite_rank = page.rank(network, damping = 0.5, personalized = initial_probabilities)$vector
network = set_vertex_attr(network, "cite_rank", value = cite_rank)
network = set_vertex_attr(network, "year", value = as.numeric(years))
whole_network = network
#remove nodes with less than 8 citations
keep = which(degree(whole_network, mode = "in") > 7) 
length(keep)
network = induced.subgraph(whole_network, keep)
top_community = which(components(network)$csize == max(components(network)$csize))
network = induced_subgraph(network, which(components(network)$membership == top_community))
attrs = get.vertex.attribute(network)
paper_ids = attrs$paper_ids
cite_rank = attrs$cite_rank
cite_rank = cite_rank/sum(cite_rank)
network = set_vertex_attr(network, "cite_rank", value = cite_rank)

#Determine the topic for each paper 
write_paper_words(field, paper_ids, "AbstractsTitleMin8")
### Run TOPIC MAPPING outside of R environment.###
topic_mapping_gammas_final = "lda_gammas_final.txt"
topic_mapping_lda_summary = "lda_summary_final.txt"
topic_mapping_word_count = "word_wn_count.txt"
num_topics = 75
topic_table = get_tt(topic_mapping_gammas_final, topic_mapping_word_count, num_topics, 10)
paper_topics = read_highest_topic(topic_mapping_gammas_final, length(paper_ids), num_topics)
network = set_vertex_attr(network, "topic", value = paper_topics)

#Determine research communities by topic and by infomap then compare them with NMI
topic_comms = create.communities(network, membership = paper_topics)
#Q = 0.504
info_comms = infomap.community(network, nb.trials = 1)
#Q = 0.474
#Get the NMI for each clustering method
nodes = vcount(network)
IM = get.vertex.attribute(network)$infomap_comm
IM = as.data.frame(cbind(1:nodes, IM))
TM = get.vertex.attribute(network)$topic
TM = as.data.frame(cbind(1:nodes, TM))
IM_TM = NMI(IM, TM)
rewired_q = get_rewired_modularity(network, 1, 50000)
#Q = 0.4257706 

network = set_vertex_attr(network, "info_comm", value = membership(info_comms))
network = set_vertex_attr(network, "cite_rank_mean", value = cite_rank) #for later manipulation
#Export the extended network
write_graph(network, "network.graphml", format = "graphml")

#Create a network where the communities are nodes
info = membership(info_comms)
communities_network = contract.vertices(network, info, vertex.attr.comb = list(paper_ids = "ignore", year = "mean", cite_rank = "sum", cite_rank_mean = gm_mean, topic = Mode, community = "ignore"))
communities_network = simplify(communities_network, remove.multiple = FALSE)
#Pretty computationally expensive, so buyer beware 
communities_network = set_edge_attr(communities_network, "weight", value = rep(1, ecount(communities_network)))
communities_network = simplify(communities_network, edge.attr.comb = list(weight = "sum"))
communities_network = set_vertex_attr(communities_network, "num_papers", value = as.numeric(sizes(info_comms)))
communities_network = set_vertex_attr(communities_network, "color", value = distinctColorPalette(vcount(communities_network)))
communities_network = set_vertex_attr(communities_network, "Topic_Str", value = sapply(get.vertex.attribute(communities_network)$topic, toString))
edge_significance = get_edge_significance(communities_network)
communities_network = set_edge_attr(communities_network, "edge significance", value = edge_significance)
write_graph(communities_network, "communities_network.graphml", format = "graphml")



## FUNCTIONS ##
get_most_connections <- function(network, TOI, p = 0.001) {
  #return the most connections between topics and topics of interest (TOI)
  #measure by sum_i^lenght(TOI) CR_TOI_i
  edge_list = get.edgelist(network)
  edge_weights = get.edge.attribute(network)$weight
  edge_significance = get.edge.attribute(network)$'edge significance'
  topics = get.vertex.attribute(network)$topic
  NOI = which(topics == TOI) #nodes of interest
  CiteRanks = get.vertex.attribute(network)$cite_rank
  num_topics = length(unique(topics))
  Topic_Sum_Connections = data.frame(cbind(0:(num_topics-1), rep(0,num_topics)))
  colnames(Topic_Sum_Connections)<-c("Topic", "Sum_Connections")
  for(i in 1:length(NOI)){
    if(i %% 10 == 0){
      print(i/length(NOI))
    }
    node = NOI[i]
    edges = which(edge_list[,1] == node)
    significant_edge_weights = sum(edge_weights[edges[which(edge_significance[edges] < p)]])
    for(edge in edges){
      significance = edge_significance[edge]
      if(significance < p){
        weight = edge_weights[edge]
        cites = edge_list[edge,2]
        topic = topics[cites]
        normal_weight = weight/significant_edge_weights
        sum_connection = normal_weight * CiteRanks[cites]
        current_sum_conn = Topic_Sum_Connections[which(Topic_Sum_Connections$Topic == topic),]$Sum_Connections
        new_sum_conn = current_sum_conn + sum_connection
        Topic_Sum_Connections[which(Topic_Sum_Connections$
                                Topic == topic),]$Sum_Connections = new_sum_conn
      }
    }
  }
  return(Topic_Sum_Connections)
}

get_edge_significance <- function(g) {
  #take a network and return a p value for each edge's significance
  heads = as.numeric(head_of(g, E(g)))
  tailz = as.numeric(tail_of(g, E(g)))
  weights = E(g)$weight
  sum_citations = 1:max(tailz)
  #get p(citing each topic)
  for(node in 1:max(tailz)){
    node_weights = weights[which(tailz == node)]
    node_weight = sum(node_weights)
    sum_citations[node] = node_weight
  }
  total_citations = sum(sum_citations)
  #get probability for each edge
  p_values = weights
  for(node in 1:max(heads)){
    edges = which(heads == node)
    node_weights = weights[edges]
    node_weight = sum(node_weights)
    for(edge in edges){
      weight = weights[edge]
      sum_cites = sum_citations[tailz[edge]]
      p_values[edge] = phyper(weight, sum_cites, total_citations - sum_cites, node_weight, lower.tail = FALSE)
    }
  }
  return(p_values)
}

read_highest_topic <- function(filename, len, num_tops){
  con = file(filename, open = "r")
  topic = vector("list", len)
  for(paper in 1:len){
    if(paper == 1000 || paper == 10000 || paper == 30000 || paper == 40000){
      cat(paper/len)
    }
    str = readLines(con, n= 1)
    lst = strsplit(str, " ")[[1]]
    lst = as.numeric(lst)
    top = which(lst == max(lst)) - 1
    topic[paper] = top
  }
  topic = as.numeric(topic)
  close(con)
  return(topic)
}

get_paper_years <- function(papers) {
  # return the year of each paper
  paper_list = papers
  num_papers = length(paper_list)
  years = dbGetQuery(db, "Select year from papers")$year
  ids = dbGetQuery(db, "Select id from papers")$id
  relevant = which(ids %in% papers)
  r_ids = ids[relevant]
  r_years = years[relevant]
  ordr = order(papers)[order(r_ids)]
  r_years = r_years[ordr]
  return(r_years)
}

get_paper_topics <- function(m){
  num_paps = length(m[,1])
  paper_topics = vector("list", num_paps)
  for (i in 1:num_paps){
    paper_topics[i] = which(m[i,] == 1) - 1
  }
  return(as.numeric(paper_topics))
}

get_true_sum_cite_ranks <- function(m, citer, num_topics){
  sum_citers = vector("list", num_topics)
  for(j in 1:num_topics) {
    sum_citers[j] = sum(m[,j]*citer)
  }
  return(as.numeric(sum_citers))
}

get_cite_sum <- function(cite_rank, comms){
  num_comms = max(comms)
  sum_cite_ranks = vector("list", num_comms)
  for(i in 1:num_comms){
    cite_rs = cite_rank[which(comms == i)]
    sum_cite_ranks[i] = sum(cite_rs)
  }
  return(as.numeric(sum_cite_ranks))
}

get_initial_probs <- function(years){
  initial_probs = pblapply(years, initial_probs_helper)
  return(as.numeric(initial_probs))
}

initial_probs_helper <- function(year) {
  exponent = (-(2017 - year))/2 #PARAMS
  return(exp(exponent))
}

get_topic_table_order <- function(top_topic, page_rnk){
  return = vector("list", max(top_topic))
  for(i in 0:max(top_topic)){
    pagers = page_rnk[which(top_topic == i)]
    return[i+1] = sum(pagers)
  }
  return(as.numeric(return))
}
print_table_latex <- function(table, lab, cap){
  xtab <- xtable(as.matrix(table), label = lab, 
                 caption = cap)
  add.to.row <- list(pos = list(0), command = NULL)
  command <- paste0("\\hline\n\\endhead\n",
                    "\\hline\n",
                    "\\multicolumn{", dim(as.matrix(table))[2] + 1, "}{l}",
                    "{\\footnotesize Continued on next page}\n",
                    "\\endfoot\n",
                    "\\endlastfoot\n")
  add.to.row$command <-command
  print(xtab, tabular.environment = 'longtable', hline.after = c(-1), add.to.row = add.to.row,
        caption.placement = "bottom")
}

get_unique_words <- function(d, n, filename){
  words_list = c()
  for(i in 1:length(d[,1])) {
    cat(i)
    words = d[i,2]
    words = strsplit(words, " ")[[1]]
    words_list = c(words, words_list)
  }
  f = file(filename, open = "r")
  word_id_count = readLines(f)
  print(length(word_id_count))
  close(f)
  first = TRUE
  for(i in 1:length(word_id_count)){
    if(i %% 5000 == 0){
      print(i/length(word_id_count))
    }
    word_id_count_i = strsplit(word_id_count[[i]], " ")[[1]]
    word = word_id_count_i[[1]]
    if(word %in% words_list){
      count = as.numeric(word_id_count_i[[3]])
      word_count = list(word, count)
      if(first == TRUE){
        word_count_d = data.frame(word_count, stringsAsFactors = FALSE)
        first = FALSE
      }
      else{
        word_count_d = rbind(word_count_d, word_count)
      }
    }
  }
  for(i in 1:length(d[,1])) {
    cat(i)
    words = d[i,2]
    words = strsplit(words, " ")[[1]]
    rank = rep(2, 100)
    for (j in 1:length(words)){
      word = words[j]
      count = word_count_d[which(word_count_d[,1] == word),2]
      rank[j] = j*as.numeric(count)
    }
    words = words[order(rank, decreasing = FALSE)]
    words = paste(words[1:n], collapse = ", ")
    d[i,2] = words
  }
  return(d)
}

get_tt <- function(filename, wn_counts, num_tops, n){
  #Create a topic table for the alzheimer's research network
  #Taking the top n unique words from filename
  f = file(filename, open = "r")
  d = data.frame(matrix(nrow = num_tops, ncol = 3), stringsAsFactors = FALSE)
  colnames(d) <- c("Topic", "Top Unique Words", "Num. of Words")
  for (i in 1:num_tops) {
    cat(i)
    for(words in 0:1){
      if(words == 0){
        title = readLines(f, n= 1)
        lst = strsplit(title, " ")[[1]]
        d[i,1] = lst[2]
        d[i,3] = lst[4]
      }
      if(words ==1){
        words = readLines(f, n = 1)
        d[i,2] = words
      }
    }
  }
  close(f)
  d = get_unique_words(d, n, wn_counts)
  return(d)
}

plot_consensus <- function(top, binary, comms){
  num_comms = max(comms)
  num_tops = length(binary[1,])
  consensus = Matrix(0, nrow = 1, ncol = num_comms)
  num_paps = consensus
  top_topic = consensus
  for(i in 1:num_comms){
    nodes = which(comms == i)
    paper_tops = binary[nodes,]
    top_count = Matrix(0, nrow = 1, ncol = num_tops)
    if(length(paper_tops) > 92) {
      for (j in 1:num_tops){
        top_count[j] = sum(paper_tops[,j])
      }
    }
    else{
      top_count = paper_tops
    }
    top_count = as.numeric(top_count)
    top_t = which(top_count == max(top_count))
    top_t = top_t[1] - 1
    ordr = order(top_count, decreasing = TRUE)
    top_count = top_count[ordr]
    paps = sum(top_count)
    con = sum(top_count[1:top])/paps
    consensus[i] = con
    num_paps[i] = paps
    top_topic[i] = top_t
  }
  consensus = as.numeric(consensus)
  num_paps = as.numeric(num_paps)
  top_topic = as.numeric(top_topic)
  hist(top_topic, "FD")
  plot(num_paps, consensus)
  return(list(num_paps, consensus, top_topic))
}

doc_assignment_binary <- function(m, len) {
  mb = Matrix(0, nrow = len, num_tops)
  for(i in 1:len){
    if(i == 100 || i == 1000 || i == 3000 || i == 7000 || i == 10000 || i == 12000){
      cat(i/len)
    }
    keepers = which(m[i,] == max(m[i,]))[1]
    Os = rep(0, num_tops)
    Os[keepers] = 1
    mb[i,] = Os
  }
  return(mb)
}

read_doc_assignment <- function(filename, len, num_tops){
  con = file(filename, open = "r")
  m = Matrix(0, nrow = len, num_tops) 
  for(paper in 1:len){
    if(paper == 1000 || paper == 10000 || paper == 30000 || paper == 40000){
      cat(paper/len)
    }
    str = readLines(con, n= 1)
    lst = strsplit(str, " ")
    lst = lapply(lst, as.numeric)[[1]]
    lst = lst/sum(lst)
    m[paper,]= lst
  }
  close(con)
  return(m)
}

create_graph_clusters <-function(edge_csv, clusters_csv) {
  # return a list that has a graph as the first element, a communities
  # object as the second and edge_table as the third
  edge_table = read.csv(edge_csv, header = FALSE)
  citation_g = graph_from_edgelist(as.matrix(edge_table))
  clusters_table = read.csv(clusters_csv, header = FALSE)
  clusters_vector = clusters_table[,2]
  # some initial readings of the communities
  # communities = make_clusters(citation_g, clusters_vector, algorithm = "SpeakEasy")
  #quick analysis
  # hist(log10(sizes(communities)), breaks = "FD", xlab = "log(community size)", main = "Community Size Distribution on Log Scale")
  return_list <- list("graph" = citation_g, "clust_t" = clusters_table, "edge_t" = edge_table)
  return(return_list)
}

community_words <- function(field, paper_list, n){
  #take a list of papers and the field we want to search, stem the words
  #return a data_frame of the top stemmed words
  words = vector("list", length(paper_list))
  for (i in 1:length(paper_list)) {
    paper = as.String(paper_list[i])
    query <- paste("SELECT", field, "FROM papers WHERE id = ", paper) 
    text <- dbGetQuery(db, query)
    corpus = word_stemmer(text)
    d <- corpus_as_data(corpus)
    #avoid double counting
    d$freq <- 1
    words[[i]] = d
  }
  bigger_d = rbind.fill(words)
  bigger_d = ddply(bigger_d, .(word), colwise(sum))
  bigger_d = arrange(bigger_d, desc(freq))
  top_d = bigger_d[which(bigger_d$freq > n),]
  return(top_d)
}

word_stemmer <- function(text) {
  #take text and return a corpus of the stemmed text and a pre-stemmed corpus
  corpus <- Corpus(VectorSource(text))
  corpus <- tm_map(corpus, content_transformer(tolower))
  corpus <- tm_map(corpus, removeWords, stopwords("english"))
  corpus <- tm_map(corpus, removeWords, c("na"))
  corpus <- tm_map(corpus, content_transformer(function(x) gsub(x, pattern = ",", replacement = " ")))
  corpus <- tm_map(corpus, removePunctuation)
  corpus <- tm_map(corpus, stripWhitespace)
  # Deal with the beta amyloid debate...
  corpus <- tm_map(corpus, content_transformer(function(x) gsub(x, pattern = "betaamyloid", replacement = "aβ")))
  corpus <- tm_map(corpus, content_transformer(function(x) gsub(x, pattern = "amyloidβ", replacement = "aβ")))
  corpus <- tm_map(corpus, content_transformer(function(x) gsub(x, pattern = "βamyloid", replacement = "aβ")))
  corpus <- tm_map(corpus, content_transformer(function(x) gsub(x, pattern = "amyloidbeta", replacement = "aβ")))
  corpus <- tm_map(corpus, content_transformer(function(x) gsub(x, pattern = "amyloid", replacement = "aβ")))
  corpus <- tm_map(corpus, content_transformer(function(x) gsub(x, pattern = "proteinaβ", replacement = "aβ")))
  corpus <- tm_map(corpus, content_transformer(function(x) gsub(x, pattern = "mouse", replacement = "mice")))
  corpus <- tm_map(corpus, content_transformer(function(x) gsub(x, pattern = "sigma1", replacement = "sigma")))
  corpus <- tm_map(corpus, content_transformer(function(x) gsub(x, pattern = "σ1", replacement = "sigma")))
  corpus <- tm_map(corpus, content_transformer(function(x) gsub(x, pattern = "mitochondria", replacement = "mitochondri")))
  corpus <- tm_map(corpus, content_transformer(function(x) gsub(x, pattern = "autophagy", replacement = "autophag")))
  corpus <- tm_map(corpus, stemDocument, language = "english")
  corpus <- tm_map(corpus, stripWhitespace)
  return(corpus)
}

write_paper_words <- function(field, papers, f) {
  # df_id_new = dataframe with id of papers and new community number
  # Return the top words in field of papers, only taking words that appear more than n times
  paper_list = papers
  f = file(f)
  lines = pblapply(paper_list, stemmed_words, field = field)
  print("writing lines...")
  lines = as.character(lines)
  writeLines(lines, f)
  close(f)
}

stemmed_words <- function(paper, field){
  #take a paper and the field we want to search, stem the words
  #return a string of the stemmed words
  query <- paste("SELECT", field, "FROM papers WHERE id = ", paper) 
  text <- dbGetQuery(db, query)
  corpus = word_stemmer(text)
  d = lapply(corpus, as.character)
  d = paste(d, collapse = " ")
  #get the year of the paper too
  return(d)
}

directed_modularity <- function(g, comms){
  #calculate the directed modularity fora given graph, g
  Q_dir = 0
  tot_edges = ecount(g)
  heads = head_of(g, E(g))
  tails = tail_of(g, E(g))
  out_degree = degree(g, mode = "out")
  in_degree = degree(g, mode = "in")
  for(comm in 1:max(comms)){
    if(comm %% 100 == 0){
      cat(comm/max(comms))
    }
    comm_nodes = which(comms == comm)
    for(head in comm_nodes){
      k_j.out = out_degree[head]
      edges = which(heads == head)
      tail_nodes = tails[edges]
      for(tail in comm_nodes){
        k_i.in = in_degree[tail]
        if(tail %in% tail_nodes){
          Q_dir = Q_dir + 1 - (k_i.in*k_j.out)/tot_edges
        }
        else{
          Q_dir = Q_dir - (k_i.in*k_j.out)/tot_edges
        }
      }      
    }
  }
  Q_dir = Q_dir/tot_edges
  return(Q_dir)
}

percent <- function(x, digits = 2, format = "f", ...) {
  paste0(formatC(100 * x, format = format, digits = digits, ...))
}

graph.reverse <- function (graph) {
  if (!is.directed(graph))
    return(graph)
  e <- get.data.frame(graph, what="edges")
  ## swap "from" & "to"
  neworder <- 1:length(e)
  neworder[1:2] <- c(2,1)
  e <- e[neworder]
  names(e) <- names(e)[neworder]
  graph = graph.data.frame(e, vertices = get.data.frame(graph, what="vertices"))
  return(graph)
}

create_network <- function() {
  #Change db into igraph format
  #Be sure to run graph.reverse script
  g = graph.reverse(graph_from_data_frame(dbReadTable(db, "paper_citation")))
  return(g)
}

get_rewired_modularity <- function(network, n, iters){
  #Rewire a given network and return the average modularity for n runs
  modularities = vector("list", n)
  for(trial in 1:n){
    cat(trial)
    rewired_net = rewire(network, with = keeping_degseq(niter = iters))
    clusters = infomap.community(rewired_net, nb.trials = 1)
    modularities[trial] = modularity(clusters)
  }
  modularities = as.numeric(modularities)
  return(modularities)
}

gm_mean = function(x, na.rm=TRUE){
  #Geometric mean 
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}



