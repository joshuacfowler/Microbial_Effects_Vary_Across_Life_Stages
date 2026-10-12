library(tidyverse)
# library(readxl)


path <- "~/Dropbox/Microbial_Effects_MetaAnalysis/raw_lit_search_files_Oct9-2026/"
files <- list.files(path = path, pattern = "savedrecs*")

  
columns_short <- c("AU", "TI", "SO", "VL", "IS", "SI", "DI", "D2", "DT", "PD", "PY", "AB", 
                   "TC", "Z9",
                   "DA", "UT")
columns <- c("Authors", "Article_Title", "Source_Title", "Volume", "Issue", "Special_Issue", "DOI",  "Book_DOI", "Document",
             "Publication_Date", "Publication_Year", "Abstract",  
             "Times_Cited_WOS_core", "Times_Cited_All_Databases", 
             "Indexed_Date", "WOS_Unique_ID")

test <- read.delim(file = paste0(path, files[7]), colClasses = "character", 
                   quote = "", 
                   row.names = NULL, stringsAsFactors = FALSE) %>% select(all_of(columns_short))



compiled_list <- list()
for(i in 1:length(files)){
  compiled_list[[i]] <- read.delim(file = paste0(path, files[i]), 
                                   colClasses = "character", quote = "", row.names = NULL, stringsAsFactors = FALSE) %>% 
    select(all_of(columns_short))
  print(i)
  }



# compiled_list[[1]]
compiled_search <- bind_rows(compiled_list) 
colnames(compiled_search) <- columns

write.csv(compiled_search, paste0(path, "WOS_search.csv"))






####### Comparing the newer search results to the previous literature extraction

effects_df <- read_csv("effects_df.csv") %>% 
  filter(experiment_id != "112-1") %>% 
  filter(!is.na(sd_RII)) %>% 
  mutate(RII = case_when(RII >= (1.0 - 1e-7) ~ 1.0,
                         RII <= (-1.0 + 1e-7) ~ -1.0, TRUE ~ RII),
         RII_01 = (RII+1)/2)

extracted_studies <- unique(effects_df$study_number)

end_of_list <- tibble(study_number = c(9453,9454,9455,9456,9457,9458,9459,9460,9461,9462,9463,9464,9465,9466,9467,9468,9469,9470,9471,9472,9473,9474,9475,9476,9477,9478,9479,9480,9481,9482,9483,9484,9485,9486,9487,9488,9489,9490,9491,9492,9493,9494,9495,9496,9497,9498,9499,9500,9501,9502,9503,9504,9505,9506,9507,9508,9509,9510,9511,9512,9513,9514,9515,9516,9517,9518,9519,9520,9521,9522,9523,9524,9525,9526,9527,9528,9529,9530,9531,9532,9533,9534,9535,9536,9537,9538,9539,9540,9541,9542,9543,9544,9545,9546,9547,9548,9549,9550,9551,9552,9553,9554,9555,9556,9557,9558,9559,9560,9561,9562,9563,9564,9565,9566,9567,9568,9569,9570,9571,9572,9573,9574,9575,9576,9577,9578,9579,9580,9581,9582,9583,9584,9585,9586,9587,9588,9589,9590,9591,9592,9593,9594,9595,9596,9597,9598,9599,9600,9601,9602,9603,9604,9605,9606),
                      include = c("dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","yes","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","yes","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped","dropped"))

path <- c("~/Dropbox/Microbial_Effects_Metaanalysis/")
old_search <- read_csv(file = paste0(path,("Microbial Effects Literature Search(WOS_all_records_2024_version).csv"))) %>% dplyr::select(c("study_number", "Authors", "Article Title", "DOI", "Publication Year", "Abstract", "UT (Unique WOS ID)")) %>% 
  rename("WOS_Unique_ID" = "UT (Unique WOS ID)") %>% 
  mutate(include = case_when(study_number %in% extracted_studies ~ "yes", study_number <=200 & !study_number %in% extracted_studies ~ "dropped", TRUE ~ NA),
         previous_study_number = study_number) %>% 
  left_join(end_of_list, join_by(study_number)) %>% mutate(include = coalesce(include.x, include.y)) %>% 
  mutate(across(everything(), as.character)) %>% 
  select(include, previous_study_number, WOS_Unique_ID)










# creating a list of the new search with  previous inclusion decision labelled out

compiled_search_1 <- compiled_search %>% left_join(old_search, by = join_by(WOS_Unique_ID)) %>% 
  mutate(study_number = row_number(),
         previous_study_number = as.numeric(previous_study_number),
         drop_notes = NA) %>% 
  dplyr::select(study_number, include, drop_notes,  Article_Title, Abstract, Authors, DOI, Source_Title, Volume, Issue, Special_Issue,  Book_DOI, Document, Publication_Date, Publication_Year, Times_Cited_WOS_core, Times_Cited_All_Databases, Indexed_Date, WOS_Unique_ID, previous_study_number)
                
write_csv(compiled_search_1, paste0(path, "Microbial_Effects_lit_filter.csv"))





# list of studies that we found good data in, but are not in this list:
found <- compiled_search_1 %>% filter(include == "yes") %>% select(previous_study_number)
not_found <- extracted_studies[!(extracted_studies %in% unique(found$previous_study_number))]



ggplot(compiled_search_1)+
  geom_point(aes(x = study_number, y = previous_study_number, color = include))+
  theme_classic()



# Noting that 49 121 122 are not included in this lit search, as far as I can tell now. Troublesome, but I guess I can just keep moving and drop those
  
  
