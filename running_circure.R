# packages
library(tidyverse)
library(data.table)
library(future)
library(future.apply)
library(dplyr)

paf_input = commandArgs(trailingOnly = TRUE)[1]
output = commandArgs(trailingOnly = TRUE)[2]

paf_file = read.delim(paf_input, header = F, sep = "\t") %>%
  rename(qname = V1,
         qlen = V2,
         qstart = V3,
         qend = V4,
         strand = V5,
         tname = V6,
         tlen = V7,
         tstart = V8,
         tend = V9,
         nmatch = V10,
         alen = V11
  )

identity = 0.95
read_map = 0.99
junction.gap = 0
max_length = 1.01
chain = 1001
chain_junction = 500 
min.overlap = 1000 
junction.min.overlap = 500 
right_extension = 1

check_0 = paf_file %>%
  ungroup() %>% mutate(check_0 = "check_0") 

check_1 = check_0 %>%
  mutate(ident = nmatch / alen) %>%
  mutate(qmap = (qend - qstart) / qlen) %>%
  filter(ident >= identity) %>%
  filter(qlen >= chain) %>%
  ungroup() %>% mutate(check_1 = "check_1") 

check_2 = check_1 %>%
  mutate(qname_tname = paste(qname, tname, sep = "__")) %>%
  mutate(spanning_read = if_else((qlen > (tlen * max_length)) &
                                   (tstart == junction.gap) &
                                   ((tlen - tend) == junction.gap) &
                                   (ident >= identity), "Y", "N")) %>%
  mutate(could_be_circular_junction = ifelse((tstart == junction.gap) |
                                               ((tlen - tend) == junction.gap), "Y", "N")) %>%  
  group_by(qname_tname) %>%
  mutate(Position.list = ifelse(could_be_circular_junction == "Y",
                                list(sort(c(tstart[could_be_circular_junction == "Y"],
                                            tend[could_be_circular_junction == "Y"]))),
                                list(NA))) %>%
  
  mutate(
    circular_junction = ifelse(
      could_be_circular_junction == "Y" &
        !is.na(Position.list) &
        length(Position.list[[1]]) == 4 &
        Position.list[[1]][1] == junction.gap &
        Position.list[[1]][4] == (tlen - junction.gap),
      "Y", "N")
  ) %>%
  mutate(
    circular_junction = ifelse(
      circular_junction == "Y" & spanning_read == "Y", 
      "N", 
      circular_junction
    )
  ) %>%
  mutate(circular_junction = ifelse(n_distinct(strand[circular_junction == "Y"]) == 1 & circular_junction == "Y", "Y", "N")) %>%
  mutate(qmap = if (any(circular_junction == "Y")) {
    sum(qmap[circular_junction == "Y"])
  } else {
    qmap
  }) %>%
  ungroup() %>%
  mutate(
    qstart_new = if_else(strand == "-", qlen - qend, qstart),
    qend_new = if_else(strand == "-", qlen - qstart, qend)
  ) %>%
  mutate(
    overhang_flag = if_else(
      (
        (qstart_new > (tstart + tlen * (max_length - 1)) & 
           tstart == junction.gap & 
           qlen > (tlen * max_length)) |
          
          ((qlen - qend_new) > ((tlen - tend) + (tlen * (max_length - 1))) & 
             (tlen - tend) == junction.gap & 
             qlen > (tlen * max_length))
      ),
      "Y", "N"
    )
  ) %>%
  select(-qstart_new, -qend_new) %>%
  mutate(
    category = case_when(
      circular_junction == "Y" ~ "Junction",   # Assign "Junction" if circular_junction == "Y"
      spanning_read == "Y" ~ "Spanning",       # Assign "Spanning" if spanning_read == "Y"
      overhang_flag == "Y" ~ "Overhang",       # Assign "Overhang" if overhang_flag == "Y"
      TRUE ~ "Normal"                          # Assign "Normal" for all other cases
    )
  ) %>%
  group_by(tname) %>%
  filter(any(category == "Junction")) %>%
  ungroup() %>% mutate(check_2 = "check_2")

check_3 = check_2 %>%
  filter((qmap >= read_map) | category == "Spanning" | category == "Overhang") %>% 
  ungroup() %>% mutate(check_3 = "check_3")

check_4 = check_3 %>%
  group_by(qname_tname) %>%
  filter(!(any(qlen > tlen  & category != "Spanning"  & category != "Overhang"))) %>% 
  ungroup() %>% mutate(check_4 = "check_4")

check_5 = check_4 %>%
  group_by(qname_tname) %>%
  mutate(group_min_chain = ifelse(category == "Junction",
                                  min((tend - tstart)[category == "Junction"]),
                                  NA_real_)) %>%
  ungroup() %>%
  filter(group_min_chain >= chain_junction | is.na(group_min_chain)) %>%
  ungroup() %>% mutate(check_5 = "check_5") %>%
  arrange(desc(category == "Junction"), tstart)

counts = check_5 %>%
  group_by(tname) %>%
  mutate(
    count_spanning_read = n_distinct(qname_tname[category == "Spanning"]),
    count_overhang_read = n_distinct(qname_tname[category == "Overhang"]),
    count_junction_read = n_distinct(qname_tname[category == "Junction"])
  ) %>%
  distinct(tname, .keep_all = T) %>%
  select(tname, count_spanning_read, count_overhang_read, count_junction_read) %>%
  ungroup()

distinct_check_0 = check_0 %>%
  distinct(tname, .keep_all = T) %>%
  select(tname, check_0)

distinct_check_1 = check_1 %>%
  distinct(tname, .keep_all = T) %>%
  select(tname, check_1)

distinct_check_2 = check_2 %>%
  distinct(tname, .keep_all = T) %>%
  select(tname, check_2)

distinct_check_3 = check_3 %>%
  group_by(tname) %>%
  filter(any(category == "Junction")) %>%
  ungroup() %>% distinct(tname, .keep_all = T) %>%
  select(tname, check_3)

distinct_check_4 = check_4 %>%
  filter((qmap >= read_map) | category == "Spanning" | category == "Overhang") %>% 
  group_by(tname) %>%
  filter(any(category == "Junction")) %>%
  ungroup() %>% distinct(tname, .keep_all = T) %>%
  select(tname, check_4)

distinct_check_5 = check_5 %>%
  filter((qmap >= read_map) | category == "Spanning" | category == "Overhang") %>% 
  group_by(tname) %>%
  filter(any(category == "Junction")) %>%
  ungroup() %>%
  distinct(tname, .keep_all = T) %>%
  select(tname, check_5)

all_checks <- reduce(list(distinct_check_0, distinct_check_1, distinct_check_2, distinct_check_3, distinct_check_4, distinct_check_5), full_join, by = "tname")

align.paf = check_5 %>%
  `colnames<-`(c("Query_name", "Query_length", "Query_start", "Query_end",
                 "Strand",
                 "Target_name", "Target_length", "Target_start", "Target_end",
                 "Number_matches", "Align_length", "check_0", "id.per.alignmnt", "Alignemnt.containment.score", "check_1", "qname_tname", 
                 "spanning_read", "could_be_circular_junction", "Position.list", "circular_junction", "overhang_flag", "category", "check_2", "check_3", "check_4", "group_min_chain", "check_5")) %>%
  select(-check_0, -check_1, -check_2, -check_3, -check_4, -check_5, -Position.list, -spanning_read, 
         -could_be_circular_junction, -circular_junction, -overhang_flag, -qname_tname, -Alignemnt.containment.score, 
         -id.per.alignmnt, -Align_length, -Number_matches, -Strand, -Query_name, -Query_start, -Query_end, -Target_length)

check.not.empty.block <- function(x){
  if(dim(x)[1]==0){
    return(F)
  } else {
    return(T)
  }
}

get.junction.reads <- function(x) {
  to_return <- x %>%
    filter(category == "Junction")
  
  if (!check.not.empty.block(to_return)) {
    return(to_return)
  }
  to_return
}

chain.support <- function(x, min.overlap, junction.min.overlap, right_extension){
  x.junctions <- x %>%
    get.junction.reads() %>%
    group_by(Target_name) %>%
    arrange(desc(group_min_chain)) %>%
    ungroup()
  
  if(!check.not.empty.block(x.junctions)){
    return(F)
  }
  
  x.not.junctions_start <- x %>%
    filter(category == "Normal") %>%
    group_by(Target_name) %>%
    arrange(desc(Query_length), desc(Target_start)) %>%
    ungroup()
  
  x.not.junctions <- x %>%
    filter(category == "Normal") %>%
    group_by(Target_name) %>%
    arrange(Target_start) %>%
    ungroup()
  
  if(!check.not.empty.block(x.not.junctions)){
    return(F)
  }
  
  
  gap.start <- x.junctions$Target_end[1]
  
  for (i in seq_len(nrow(x.not.junctions_start))) {
    
    #print(paste("Checking row", i, 
    #            "x.not.junctions_start$Target_start:", x.not.junctions_start$Target_start[i], 
    #            "x.not.junctions_start$Target_end:", x.not.junctions_start$Target_end[i]))
    
    if (x.not.junctions_start$Target_start[i] <= x.junctions$Target_end[1] - junction.min.overlap &&
        x.not.junctions_start$Target_end[i] > x.junctions$Target_end[1] + junction.min.overlap + right_extension) {
      
      gap.start <- x.not.junctions_start$Target_end[i]
      #print(paste("Match found at row", i, "Updated gap.start:", gap.start))
      break 
    }
  }
  
  #print(paste("Final gap.start:", gap.start))
  
  gap.end <- x.junctions$Target_start[2] 
  
  
  #print(x.junctions[1:2,])
  #print(c(gap.start, gap.end))
  
  while(((gap.start - gap.end) < min.overlap) & check.not.empty.block(x.not.junctions)){
    #print(dim(x.not.junctions))
    #print(((gap.start - gap.end) < min.overlap) || !check.not.empty.block(x.not.junctions))
    if((gap.start - x.not.junctions$Target_start[1]) >= min.overlap){
      if((x.not.junctions$Target_end[1] - gap.start) >= right_extension){
        gap.start <- x.not.junctions$Target_end[1]
        #print(x.not.junctions[1,])
      }
    }
    x.not.junctions <- x.not.junctions[-1,]
    #print(c(gap.start, gap.end))
  }
  
  if((gap.start - gap.end) > junction.min.overlap){
    return(T)
  } else {
    return(F)
  }
}


is.contig.chained <- function(x, min.overlap, junction.min.overlap, right_extension){
  #print(x$Target_name[1])
  
  to_return <- chain.support(x, min.overlap, junction.min.overlap, right_extension)
  return(to_return)
}

options(future.globals.maxSize = 8000 * 1024^2)
plan(multicore)

prediction <- future_lapply(unique(align.paf$Target_name), function(x) {
  is.contig.chained(align.paf %>% filter(Target_name == x), min.overlap, junction.min.overlap, right_extension)
})

plan(sequential)

chained.contigs <- data.frame(
  tname = unique(align.paf$Target_name), 
  prediction = unlist(prediction)
) %>%
  #filter(prediction == T) %>%
  right_join(all_checks, by = "tname") %>%
  left_join(counts, by = "tname") %>%
  mutate(sample = paf_input) %>%
  mutate(sample = sub(".*/(.*)\\.paf$", "\\1", sample)) %>%
  select(tname, prediction, count_spanning_read, count_overhang_read, count_junction_read, sample) %>%
  mutate(prediction = if_else(is.na(prediction), FALSE, prediction)) %>%
  #filter(prediction == TRUE) %>%
  rename(contig = tname,
          reads_longer_than_contig_no_ab_split = count_spanning_read,
          reads_overhanging = count_overhang_read,
          reads_mapping_over_ab = count_junction_read,
          file = sample) %>%
  select(contig, file, prediction, reads_mapping_over_ab, reads_longer_than_contig_no_ab_split, reads_overhanging)
  

write.table(chained.contigs, output, sep="\t", row.names=F, col.names=T, quote=F)