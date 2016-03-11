# this is an R script
# Create a block of text to use for uploading custom tracks to UCSC; this version is customized for bigWig files with hg19 genome


# Upload the files here:
# https://genome.ucsc.edu/cgi-bin/hgCustom
# more instructions here: # http://genome.ucsc.edu/goldenPath/help/bigWig.html

# sample : 
#  track type=bigWig name="My Big Wig" description="A Graph of Data from My Lab" bigDataUrl=http://myorg.edu/mylab/myBigWig.bw 

# Custom function
# function to paste the information together into a block of text that we can copy/paste into UCSC to upload a custom track

UCSC_track_string_func<-function(file_name,URL_base,Track_Type){
  file_URL<-paste0(URL_base,file_name)
#   x<-paste0("track type=",Track_Type," name=",file_name," dataUrl=",file_URL) # for a wig
    x<-paste0("track type=",Track_Type," name=",file_name," bigDataUrl=",file_URL) # for a bigwig

  return(x)
}


# the file path to the dir containing the bigwigs; accessible via http
Wig_Dir<-"/data/sequence/bigwigs/"

# get a list of all the files in the dir; dir should contain ONLY the files to be uploaded to UCSC, all ofthe same type, bigWig in this case
Wig_List<-dir(Wig_Dir)

# sample Wig_List entries look like this:
# "ABC-F-H3K27AC.bw" ; H3K27AC is a histone mark that we care about, ABC is the sample name, F is the treatment

strsplit(x = Wig_List[1],split = "-")[[1]][1] # sample
strsplit(x = Wig_List[1],split = "-")[[1]][3] # hist mark
# H3K27AC
# H3K27ME3
# H3K4ME3
# H3K9AC
# H3K9ME3
# INPUT

# make a df with a col for the hist mark
Wig_df<-data.frame(File=Wig_List,
                   Hist=gsub(pattern = ".bw",replacement = "",x =sapply(X = Wig_List,FUN = function(x)  strsplit(x = x,split = "-")[[1]][3]))
                   )


# we want to reorder them all by the histone mark so they show up in order in UCSC
Wig_df<-Wig_df[with(Wig_df, order(Hist)), ]

# save a vector of entries
Wig_List_sortHist<-as.character(Wig_df[["File"]])
class(Wig_List_sortHist)




# the URL to the files on the server
Wig_URL_base<-"https://internets.com/sequence/bigwigs/"
Wig_URL_list<-paste0(rep(Wig_URL_base,times=length(Wig_List_sortHist)),Wig_List_sortHist)
Track_Type<-"bigWig"

# preview the string to be created
UCSC_track_string_func(Wig_List_sortHist[1],Wig_URL_base,Track_Type)

# apply the function for all of the files
UCSC_track_string<-sapply(Wig_List_sortHist,UCSC_track_string_func,URL_base=Wig_URL_base,Track_Type=Track_Type)
UCSC_entry<-paste(UCSC_track_string,Wig_URL_list,sep = "\n")

# this prints text to the console in the format we need to copy/paste into UCSC
cat(UCSC_entry)





# make sure to save this in a Session:
# https://genome.ucsc.edu/cgi-bin/hgSession
# http://genome.ucsc.edu/goldenPath/help/bigWig.html
