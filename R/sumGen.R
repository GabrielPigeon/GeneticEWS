require(diveRsity)


genopop_81_10 <- read.csv (file = "~/doctorat/RAnalyse/GenEWS/data/Genepop_RAM_1981_2010.csv", header = FALSE)
basic_ram <- divBasic(genopop_81_10)
RamGen <- data.frame(yr=colnames(basic_ram$Allelic_richness),
  # He=basic_ram$He["overall",],
                     Ar=basic_ram$Allelic_richness["overall",]
)
basic_ram$Allele_number
basic_ram$locus_pop_size  
# basicStats(genopop_81_10)
  
basic_ram$Allelic_richness
rg2 <- readxl::read_xlsx("~/doctorat/RAnalyse/GenEWS/data/summary_1981_2009_GenAlex.xlsx")
rg2$YEAR
RamGen <- merge(RamGen,rg2[,c("YEAR","He","Fis")],by.x = "yr",by.y = "YEAR")
write_csv(RamGen,"../EarlySignal/data/ram/genetics.csv")
