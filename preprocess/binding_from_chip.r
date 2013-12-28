# master_dat is a table of TF names and the corresponding file names (since ENCODE is lovely and these are all nigh-on incomprehensible), assuming these files are in BigBed format...
master_dat<-read.table('tf_names_and_files.txt',as.is=TRUE)
# yes, this is just the NAME of the list
# ideally this is the DEFINITIVE list...
peak_list<-'peak_list'

N_FACTORS<-nrow(master_dat)
for (i in 1:N_FACTORS){
    tf_name<-master_dat[i,1]
    tf_filename<-master_dat[i,2]
    print(tf_name)
    print(tf_filename)

    make_bedGraph<-paste("./bigBedToBed ",tf_filename, " temp_chip.bed",sep="")
    system(make_bedGraph)
    cat("Obtaining bound information!\n")
    get_bound_peaks<- paste("bedmap --indicator --delim '\t' ",peak_list," temp_chip.bed > boundpeaks.bed",sep="")
    system(get_bound_peaks)
    cat("Commencing data processing!\n")

    add_name <- paste("awk 'BEGIN{ print \"",tf_name,"\" }{ print }' boundpeaks.bed > boundpeaks_withname.bed",sep="")
    system(add_name)
    if (i==1){
        system("mv boundpeaks_withname.bed binding_mat")
    }
    else{
        forming_matrix <- "paste binding_mat boundpeaks_withname.bed > temp_mat"
        system(forming_matrix)
    }
    rename_temps <- "mv temp_mat binding_mat"
    system(rename_temps)

    cleanup<-"rm temp_chip.bed boundpeaks.bed boundpeaks_withname.bed"
    system(cleanup)
}
