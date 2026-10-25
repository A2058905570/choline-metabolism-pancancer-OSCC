rm(list = ls())
setwd('/data/nas3/mayanyu/project/06_TJZK-10306-1/R/')

path <- getwd()
list <- list.files(path = './')
for (i in list) {

  script_path <- i
  script_content <- readLines(script_path, warn = FALSE)
  

  remove_comments <- function(script_lines) {

    script_lines <- gsub("^\\s*
    return(script_lines)
  }
  

  script_without_comments <- remove_comments(script_content)
  

  writeLines(script_without_comments, script_path)
  
  
}
rm(list = ls())
setwd('/data/nas3/mayanyu/project/06_TJZK-10306-1/R/')




path <- getwd()
list <- list.files(path = './')
i <- 'r.04_DEG_tumor.r'
for (i in list) {

  script_path <- i
  script_content <- readLines(script_path, warn = FALSE)
  
  j <- 30
  for (j in 1:length(script_content)) {

    first_line <- script_content[j]
    
    num_hashtags <- sum(gregexpr("#(?![#])", first_line, perl = TRUE)[[1]] >= 0)
    if(num_hashtags > 1){next}
    

    
    clean_code <- gsub("



    print(clean_code)
    script_content[j] <- clean_code
    writeLines(script_content, script_path)
    
  }
  
  
}
