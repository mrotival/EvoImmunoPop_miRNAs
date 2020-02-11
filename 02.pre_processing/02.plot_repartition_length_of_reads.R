#########################################################################
##This script aims to plots severals informations such as the number of##
##read per samples, repartition of length of the samples etc ...       ##
#########################################################################
#############
##Libraries##
#############
require(data.table)

#############
##load data##
#############
samples_informations = fread(paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/01.description_of_raw_data/locations_and_informations_of_fastq_files.tsv", sep=""))

###################
##Number of Reads##
###################
get_number_of_reads <- function(individual, condition, library_number, batch_name){
  temp = fread(paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/02.pre_processing/number_different_reads/", batch_name, "_", individual,"-", condition, "_",  library_number, ".all_reads.count.tsv", sep=""))
  temp[, individual := individual]
  temp[, condition := condition]
  temp[, library_number := library_number]
  temp[, batch_name := batch_name]
  return(temp)
}

all_reads_data = list()
for (i in samples_informations[, c(1:.N)]){
  all_reads_data[[i]] = get_number_of_reads(samples_informations[i, individual], samples_informations[i, condition], samples_informations[i, library_number], samples_informations[i, batch_name])
}

all_reads_data = rbindlist(all_reads_data)
print(all_reads_data)


##################################
##plot number of reads per fastq##
##################################
to_plot = all_reads_data[1:.N]
to_plot[, sample := paste(individual, condition, library_number, batch_name, sep="_")]
to_plot  = to_plot[, .(number_of_read = sum(number_of_read)), by = sample]
to_plot = to_plot[ order(number_of_read)]
to_plot[, sample := factor(sample, levels = unique(sample))]

require(ggplot2)
p <- ggplot(to_plot, aes(x = sample, y = number_of_read, color = NULL))
p <- p + geom_bar(stat = "identity", fill = "Blue", alpha = 0.5)
p <- p + theme_classic()
p <- p + theme(axis.text.x=element_blank(),
                axis.ticks.x=element_blank())
p <- p + xlab("fatsq")
p <- p + ylab("number of read in fastq file")

pdf(paste(EVO_IMMUNO_POP,"Maxime/miRNA_V2/figures/02.pre_processing/number_of_reads_per_fastq.pdf",sep=''), width = 5, height = 5)
print(p)
dev.off()

###################################
##plot number of reads per sample##
###################################
to_plot = all_reads_data[1:.N]
to_plot[, sample := paste(individual, condition, library_number,sep="_")]
to_plot  = to_plot[, .(number_of_read = sum(number_of_read)), by = sample]
to_plot = to_plot[ order(number_of_read)]
to_plot[, sample := factor(sample, levels = unique(sample))]

require(ggplot2)
p <- ggplot(to_plot, aes(x = sample, y = number_of_read, color = NULL))
p <- p + geom_bar(stat = "identity", fill = "Blue", alpha = 0.3)
p <- p + theme_classic()
p <- p + theme(axis.text.x=element_blank(),
                axis.ticks.x=element_blank())
p <- p + xlab("sample")
p <- p + ylab("number of read in fastq file")

pdf(paste(EVO_IMMUNO_POP,"Maxime/miRNA_V2/figures/02.pre_processing/number_of_reads_per_sample.pdf",sep=""), width = 5, height = 5)
print(p)
dev.off()

print(to_plot[, sum(number_of_read>8e6)/.N])

##################################
##get quality of reads per fastq##
##################################
get_quality_of_reads <- function(individual, condition, library_number, batch_name){
  temp = fread(paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/02.pre_processing/QC/", batch_name, "_", individual,"-", condition, "_",  library_number, ".quality.txt", sep=""))
  temp[, individual := individual]
  temp[, condition := condition]
  temp[, library_number := library_number]
  temp[, batch_name := batch_name]
  temp[, sum := as.numeric(sum)]
  return(temp)
}

quality_data = list()
for (i in samples_informations[, c(1:.N)]){
  quality_data[[i]] = get_quality_of_reads(samples_informations[i, individual], samples_informations[i, condition], samples_informations[i, library_number], samples_informations[i, batch_name])
}
quality_data = rbindlist(quality_data)


######################
##plot QC per sample##
######################
to_plot = quality_data[1:.N]
to_plot[, sample := paste(individual, condition, library_number, batch_name,sep="_")]
samples = to_plot[, unique(sample)]

# for (i in 1:110){
#   print(i)
#   bloup = to_plot[ sample %in% samples[1+(10*(i-1)):min(10+(10*(i-1)), length(samples))]]
#   # print(bloup[, length(unique(sample))])
#   p <- ggplot(bloup, aes(x = column, y = mean, color = sample))
#   p <- p + geom_line() +geom_point()
#   p <- p + theme_bw()
#   pdf(paste(EVO_IMMUNO_POP,"Maxime/miRNA_V2/figures/02.pre_processing/qcplots/qcplots_", i, ".qcplots.pdf", sep=""), width = 8, height = 3)
#   print(p)
#   dev.off()
# }


##case of AFB040-1
bloup = to_plot[ (individual == "AFB040") & (condition == 1)]

p <- ggplot(bloup, aes(x = column, y = mean, color = sample))
p <- p + geom_line() +geom_point()
p <- p + theme_bw()
pdf(paste(EVO_IMMUNO_POP,"Maxime/miRNA_V2/figures/02.pre_processing/AFB040.pdf", sep=""), width = 8, height = 3)
print(p)
dev.off()


#########################
##Length after clipping##
#########################
get_number_of_reads_clipped <- function(individual, condition, library_number, batch_name){
  temp = fread(paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/02.pre_processing/number_different_reads/", batch_name, "_", individual,"-", condition, "_",  library_number, ".clipped.count.tsv", sep=""))
  temp[, individual := individual]
  temp[, condition := condition]
  temp[, library_number := library_number]
  temp[, batch_name := batch_name]
  return(temp)
}

clipped_reads_data = list()
for (i in samples_informations[, c(1:.N)]){
  clipped_reads_data[[i]] = get_number_of_reads_clipped(samples_informations[i, individual], samples_informations[i, condition], samples_informations[i, library_number], samples_informations[i, batch_name])
}
clipped_reads_data = rbindlist(clipped_reads_data)


########################
##plot length of reads##
########################
to_plot = clipped_reads_data[1:.N]
to_plot[, sample := paste(individual, condition, library_number, batch_name,sep="_")]
samples = to_plot[, unique(sample)]
to_plot = to_plot[, .(length_of_read = length_of_read, frequency = number_of_read /sum(number_of_read)), by = c("sample")]
for (i in 1:110){
  print(i)
  bloup = to_plot[ sample %in% samples[1+(10*(i-1)):min(10+(10*(i-1)), length(samples))]]
  # print(bloup[, length(unique(sample))])
  p <- ggplot(bloup, aes(x = length_of_read, y = frequency, color = sample))
  p <- p + geom_line() +geom_point()
  p <- p + theme_bw()
  pdf(paste(EVO_IMMUNO_POP,"Maxime/miRNA_V2/figures/02.pre_processing/length_clipped_reads/length_of_clipped_reads_", i, ".length.pdf", sep=""), width = 8, height = 3)
  print(p)
  dev.off()
}
