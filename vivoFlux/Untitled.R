lo = d.normalized.tidy %>% 
  filter(phenotype == "WT" & infused.tracer == "Glucose" & Compound == "Glucose")

lo %>% 
  ggplot(aes(x = sample, y = enrichment, fill = as.character(C_Label))) + 
  geom_bar(stat = "identity", position = "stack")  +
  scale_y_continuous(breaks = seq(0, 1, .1)) +
  scale_fill_brewer(palette = "Dark2") +
  geom_hline(yintercept = .18)

lo %>% tail()
lo.1 =  lo %>% filter(sample == "a-art-L4")
lo.1 =  lo %>% filter(sample == "aa-tail-L8")
lo.1 =  lo %>% filter(sample == "z-tail-L7")

((lo.1$C_Label) / (lo.1$C_Label.max) * (lo.1$enrichment) ) %>% sum()

plt.enrich.molecules.glucose + geom_hline(yintercept = .9, color = "red")


d.normalized.tidy %>% 
  filter(sample == "a-art-L4") 
