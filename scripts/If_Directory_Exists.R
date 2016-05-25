# Check to see if directory exists. If so, delete the directory and re-create it. If not, make the directory

# the directories we are using
DIRS<-c("/home/user/DataOutput/",
        "/home/user/DataOutput/A/",
        "/home/user/DataOutput/B/")

# iterate through the vector of directories
for(i in 1:length(DIRS)){
  if(dir.exists(DIRS[i])){
  # delete the directory
    unlink(DIRS[i],recursive = T,force = T)
    
    # make a new copy of the directory
    dir.create(DIRS[i],mode = "0770")
  } else {dir.create(DIRS[i],mode = "0770")}
} 

# when running scripts, we need to make sure old results are removed
# or they might accidentally be used in place of new data if the script fails
# to do this, delete the data output directories we will be working with
# and make them over again each time the script is run
# this runs in R


# this works well when these commands are at the start of a script
