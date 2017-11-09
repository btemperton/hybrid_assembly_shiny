clusters <- df %>% 
  group_by(cluster_name) %>%
  summarise(representative_length=max(contig_len),
            number_of_members=n(),
            mean_pct_id=mean(similarity_to_cluster_representative, na.rm=TRUE),
            representative=contig_id[which.max(contig_len)],
            type=substr(contig_id[which.max(contig_len)], 1, 1)) %>%
  arrange(desc(representative_length))
write.csv(clusters, 'clusters.csv', row.names = FALSE, quote=FALSE)
