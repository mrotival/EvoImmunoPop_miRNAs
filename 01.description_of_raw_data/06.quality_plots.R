###########
##imports##
###########
require(ggplot2)


###############################
## Number of read per sample ##
###############################
input_data = fread(paste(EVO_IMMUNO_POP, "/Maxime/miRNA_V2/data/01.description_of_raw_data/locations_and_informations_of_fastq_files.tsv", sep=""))
input_data[, index := NULL]

get_number_of_read <- function(i){
  print(i)
  #Number Complete reads
  temp_result = input_data[i, .(individual,  condition, library_number, batch_name, library_ID)]
  complete_reads_info = fread(paste(EVO_IMMUNO_POP, "/Maxime/miRNA_V2/data/02.pre_processing/number_different_reads/", temp_result$batch_name[1],"_", temp_result$library_ID[1], ".all_reads.count.tsv", sep=""))
  names(complete_reads_info) = paste(names(complete_reads_info), "all_reads", sep="..")
  results = cbind(temp_result, complete_reads_info)
  return(results)
}

to_plot = list()
for (i in 1:nrow(input_data)){
  to_plot[[i]] = get_number_of_read(i)
}
to_plot = rbindlist(to_plot)
to_plot = to_plot[order(number_of_read..all_reads)]
to_plot[, axis := 1:.N]
p <- ggplot(to_plot, aes(x = axis, ymax = number_of_read..all_reads))
p <- p + geom_ribbon(ymin=0, fill = "#0db8c1", alpha = 0.5, color = "#0e79a3")
p <- p + theme_bw()
p <- p + theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        text = element_text(size=14))
p <- p + ylab("Number of raw reads per fastq files")

pdf(sprintf("%s/Maxime/miRNA_V2/figures/01.description_of_raw_data/Number_of_raw_reads_per_fastq.pdf",EVO_IMMUNO_POP), width = 4, height =4)
print(p)
dev.off()


to_plot = to_plot[number_of_read..all_reads<=2e6]
to_plot[, library_ID := factor(library_ID, levels = unique(library_ID))]

p <- ggplot(to_plot, aes(x = library_ID, y = number_of_read..all_reads))
p <- p + geom_point(color = "#0e79a3")
p <- p + theme_bw()
p <- p + theme(axis.title.x=element_blank(),
        text = element_text(size=12),
        axis.text.x = element_text(angle = 60, hjust = 1))
p <- p + ylab("Number of read in the lowast fastq")

pdf(sprintf("%s/Maxime/miRNA_V2/figures/01.description_of_raw_data/Number_of_raw_reads_per_fastq_problematics_samples.pdf",EVO_IMMUNO_POP), width = 4, height =4)
print(p)
dev.off()


#############################
## Number of clipped reads ##
#############################

get_number_of_clipped_read <- function(i){
  print(i)
  #Number Complete reads
  temp_result = input_data[i, .(individual,  condition, library_number, batch_name, library_ID)]
  complete_reads_info = fread(paste(EVO_IMMUNO_POP, "/Maxime/miRNA_V2/data/02.pre_processing/number_different_reads/", temp_result$batch_name[1],"_", temp_result$library_ID[1], ".clipped.count.tsv", sep=""))
  names(complete_reads_info) = paste(names(complete_reads_info), "clipped", sep="..")
  results = cbind(temp_result, complete_reads_info)
  return(results)
}

to_plot = list()
for (i in 1:nrow(input_data)){
  to_plot[[i]] = get_number_of_clipped_read(i)
}
to_plot = rbindlist(to_plot)
to_plot = to_plot[(length_of_read..clipped>=18) & (length_of_read..clipped<=26)]
to_plot = to_plot[, .(number_read_clipped_size_filtered = sum(number_of_read..clipped)), by = c("batch_name", "library_ID")]


to_plot = to_plot[order(number_read_clipped_size_filtered)]
to_plot[, axis := 1:.N]
p <- ggplot(to_plot, aes(x = axis, ymax = number_read_clipped_size_filtered))
p <- p + geom_ribbon(ymin=0, fill = "#0db8c1", alpha = 0.5, color = "#0e79a3")
p <- p + theme_bw()
p <- p + theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        text = element_text(size=14))
p <- p + ylab("Number of clipped reads of good size")

pdf(sprintf("%s/Maxime/miRNA_V2/figures/01.description_of_raw_data/Number_of_reads_clipped_filtered1826.pdf",EVO_IMMUNO_POP), width = 4, height =4)
print(p)
dev.off()


##############
##GC percent##
##############
get_GC_percent <- function(i){
  print(i)
  #Number Complete reads
  temp_result = input_data[i, .(individual,  condition, library_number, batch_name, library_ID)]
  complete_reads_info = fread(paste(EVO_IMMUNO_POP, "Maxime/miRNA_V2/data/02.pre_processing/QC/", temp_result$batch_name[1],"_", temp_result$library_ID[1], ".quality.txt", sep=""))
  gc_percent = complete_reads_info[, (sum(as.numeric(G_Count)) + sum(as.numeric(C_Count)))/sum(as.numeric(Max_count)) ]
  temp_result[, gc_percent := gc_percent]
  results = temp_result[1:.N]
  return(results)
}

to_plot = list()
for (i in 1:nrow(input_data)){
  to_plot[[i]] = get_GC_percent(i)
}
to_plot = rbindlist(to_plot)
to_plot = to_plot[order(gc_percent)]
to_plot[, axis := 1:.N]
p <- ggplot(to_plot, aes(x = axis, y = gc_percent))
p <- p + geom_line(color = "#0e79a3")
p <- p + theme_bw()
p <- p + theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        text = element_text(size=14))
p <- p + ylab("GC percent")

pdf(sprintf("%s/Maxime/miRNA_V2/figures/01.description_of_raw_data/GC_percent.pdf",EVO_IMMUNO_POP), width = 4, height =4)
print(p)
dev.off()
