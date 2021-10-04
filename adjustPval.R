#!/biosw/debian7-x86_64/R/3.2.2/bin/Rscript --vanilla

## Retrieve arguments
args=commandArgs(TRUE)

## Smallest number before p-values are set to 0
## this is necessary as other programs (e.g. awk) report WRONG results
MINNUM=2.2250738585072014e-308

## Help
help <- function(){
  cat("\nAdjust p-values\n")
  cat("Usage: adjustPval.R -i - -c N -o -\n")
  cat("-i : input table as a file or stdin (-)\n")
  cat("-c : column with raw p-values (int)\n")
  cat("-o : Output as a file or stdout (-)\n")
  cat("\n")
  q()
}

## Save values of each argument
if(length(args)==0 || !is.na(charmatch("-help",args))){
  help()
} else {
  for(ii in 1:length(args)){
    if(grepl("^-",args[ii]) && args[ii] != "-"){
      if(ii+1<=length(args) && (!grepl("^-",args[ii+1]) || args[ii+1]=="-")){
        assign(gsub("-","",args[ii]),args[ii+1])
      } else {assign(gsub("-","",args[ii]),1) }
    }
  }}


## Set the ouput path : File or STDOUT
if(exists("o")){
  if(o=="stdout" || o=="-"){
    output=stdout()
  } else {output=o}
} else { output=stdout() }

## Check if column is specified
if(exists("c")){
  col=strtoi(c)
  if(is.na(col) || col<1){ cat("P-value column is not int>=1 (-c)\n"); q() }
} else { cat("P-value column not specified (-c)\n"); q() }

## Load the matrix into a dataframe
if(exists("i")){
  if(i=="stdin" || i=="-"){
    count_table=read.csv(pipe('cat /dev/stdin'), sep="\t", skip=0, header = F, comment.char = "", check.names = F)
  } else if (file.exists(i)){
    count_table=read.csv(i, sep="\t", skip=0, header = F, comment.char = "", check.names = F)
  } else { cat("Input file does not exist\n"); q() }
  
  ##Test if the matrix has at least col columns
  if(ncol(count_table)<col){
    cat("The matrix does not have enough columns\n"); q()
  }
  
  ## Test column c to see if they contain only floats
  if(!all(sapply(count_table, function(x) class(x) %in% c("integer", "numeric"))[col])){
    cat("Column c does not contain only numbers\n");q()
  }

} else { cat("No input specified\n"); q() }


## Adjust p-value using BH
adj_pval <- p.adjust(count_table[[col]], method = "BH")
adj_pval[adj_pval<=MINNUM]=0

## append adj p-value column to table
count_table=cbind(count_table,adj_pval)

## Writing the results
write.table(count_table,
            sep="\t",
            quote=FALSE, row.names=FALSE, col.names=FALSE)

