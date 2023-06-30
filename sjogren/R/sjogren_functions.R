#SJOGREN PACKAGE FUNCTIONS

#' Clonotyping function
#'
#' This function clonotypes VDJ contig files based on identical heavy chain VDJ assignments, identical light chain VJ assignments,
#' #and a user-defined threhold for the amount of nucleotide mismatches in the CDR3 region.
#' @param contig The VDJ contigs (as a data frame, containing at least a barcode, VDJ assignments, and CDR3 nucleotide sequences).
#' @param threshold The allowed number of mismatches (as an integer) in the CDR3 region between two sequences in a clonotype.
#' @param levenshtein Should the distance be calculated based on the levenshtein distance (allowing indels)? Defaults to FALSE (using point mutation mismatches instead)
#' @keywords clonotyping
#' @returns The original contig file, with clonotype numbers for each cell in the "clonotype_new" column.
#' @export
#' @examples contig <- clonotype(contig, threshold = 2, levenshtein = F)

clonotype <- function(contig, threshold, levenshtein = F){

  #Make sure correct packages are installed
  stopifnot(require(fastmatch))
  stopifnot(require(stringi))
  stopifnot(require(stringr))
  stopifnot(require(udpipe))
  stopifnot(require(usefun))
  stopifnot(require(dplyr))
  stopifnot(require(tidyr))

  `%fin%` <- function(x, table) {
    stopifnot(require(fastmatch))
    fmatch(x, table, nomatch = 0L) > 0L
  }

  #Require a donor and cell barcode column
  if (!(all(c("donor", "barcode", "ident") %in% colnames(contig)))){
    stop("The object should have a column stating the donor and identity (ident), and it should have barcodes.")
  }

  if (!(c("cell_barcode") %in% colnames(contig))){
    contig$cell_barcode <- paste0(contig$ident, "_", contig$barcode)
  }

  #Initialisation message
  print("Initialising...")

  #Only keep productive chains
  contig_beforefilter <- contig #save the starting contig table

  contig <- contig[contig$cdr3_nt != "None" & contig$productive %in% c("True", "true", T),]

  contig_nonproductive <- contig_beforefilter[contig_beforefilter$cdr3_nt == "None" | !(contig_beforefilter$productive %in% c("True", "true", T)),]

  print(paste0(length(unique(contig_nonproductive)), " cells had chains removed due to them being unproductive." ))

  #make sure threshold is okay
  if (threshold > min(nchar(contig$cdr3_nt)) | threshold < 0 | !is.numeric(threshold)){
    stop("Please choose an appropriate threshold. This error can occur if the threshold is below 0, more than the minimum CDR3 nucleotide length, or if the threshold is non-numeric.")
  }



  contig <- data.table(contig) #convert to data.table for speed
  contig_start <- contig #save the starting contig table
  suppressWarnings(contig$clonotype_new <- NULL) #remove any pre-rexisting clonotype columns

  #If a gene is stated as "None", set to empty
  contig$d_gene[contig$d_gene == "None"] <- ""

  #Take out all IGH chains
  contig_IGH <- contig[contig$chain == "IGH",]
  #Remove cells with dual IGH chains
  dual_chains <- contig_IGH$cell_barcode[duplicated(contig_IGH$cell_barcode)]
  contig_IGH <- contig_IGH[!(contig_IGH$cell_barcode %fin% dual_chains),]

  print(paste0(length(dual_chains), " cells with two IGH chains were found."))


  #Group the cells with the same IGH genes
  contig_IGH_grouping <- contig_IGH %>%
    group_by(v_gene, d_gene, j_gene) %>%
    dplyr::count()

  clonotypes_IGH <- seq(1:nrow(contig_IGH_grouping)) #define IGH clonotypes, which is the amount of unique IGH genes

  contig_IGH_grouping$n <-NULL
  contig_IGH_grouping$clonotype_IGH <- clonotypes_IGH #add the IGH clonotypes to the genes

  contig <- left_join(contig, contig_IGH_grouping, by = c("v_gene", "d_gene", "j_gene")) #add the IGH clonotypes to the cells


  #Group cells with the same IGKL chains, same procedure as for IGH
  contig_IGKL <- contig[contig$chain == "IGK" | contig$chain == "IGL",]


  contig_IGKL_grouping <- contig_IGKL %>%
    group_by(v_gene, j_gene) %>%
    dplyr::count()

  max_IGH <- max(clonotypes_IGH)+1
  max_IGKL <- max(clonotypes_IGH) + nrow(contig_IGKL_grouping) #the IGKL clonotype numbering starts where the IGH clonotype numbering left off

  clonotypes_IGKL <- seq(max_IGH, max_IGKL)

  contig_IGKL_grouping$n <-NULL
  contig_IGKL_grouping$clonotype_IGKL <- clonotypes_IGKL

  contig <- left_join(contig, contig_IGKL_grouping, by = c("v_gene", "j_gene"))


  #set the clonotypes of dual IGH chains to "none"
  contig$clonotype_IGH[contig$cell_barcode %fin% dual_chains] <- "none, dual IGH"
  contig$clonotype_IGKL[contig$cell_barcode %fin% dual_chains] <- "none, dual IGH"

  #check for dual or missing IGH/IGKL chains
  check_IG <- contig %>%
    group_by(cell_barcode, chain) %>%
    dplyr::count()

  #also set the clonotype of cells that have only one chain (IGH or IGKL) to "none"
  check_IGH <- check_IG[check_IG$chain == "IGH",]
  check_IGKL <- check_IG[check_IG$chain == "IGK" | check_IG$chain == "IGL",]

  one_chain <- outersect(check_IGH$cell_barcode, check_IGKL$cell_barcode) #outersect looks for cells that only appear in one of the two data frames
  contig$clonotype_IGH[contig$cell_barcode %fin% one_chain] <- "none, only one chain"
  contig$clonotype_IGKL[contig$cell_barcode %fin% one_chain] <- "none, only one chain"

  print(paste0(length(one_chain), " cells with only one chain (IGH or IGK/L) were found."))

  #For dual IGK and IGL chains in one cell, keep the entry with the highest amount of reads
  contig_IGKL <- contig[contig$chain == "IGK"| contig$chain == "IGL",]
  contig_IGKL_dup <- contig_IGKL[contig_IGKL$clonotype_IGKL != "none, dual IGH" & contig_IGKL$clonotype_IGKL != "none, only one chain",]

  contig_IGKL_dup <- contig_IGKL_dup[duplicated(contig_IGKL_dup$cell_barcode),] #cells with two IGKL chains will appear twice in the data frame
  contig_IGKL_dup <- contig_IGKL_dup %>%
    group_by(cell_barcode) %>% top_n(1, reads) #keep the entry with the highest reads

  if (nrow(contig_IGKL_dup) > 0){
    contig_IGH_counterparts <- contig[contig$cell_barcode %fin% contig_IGKL_dup$cell_barcode,]
    contig_IGH_counterparts <- contig_IGH_counterparts[contig_IGH_counterparts$chain == "IGH",]
    contig_IGKL_dup <- rbind(contig_IGKL_dup, contig_IGH_counterparts)
    contig <- contig[!(contig$cell_barcode %fin% contig_IGKL_dup$cell_barcode),] #replace the dual IGKL entries with the one with the most reads
    contig <- rbind(contig, contig_IGKL_dup)
  }


  #Group the cells based on their IGKL and IGH chains
  contig_full_grouping <- contig %>%
    group_by(cell_barcode, clonotype_IGH, clonotype_IGKL) %>%
    dplyr::count()

  contig_full_grouping$clonotype_IGH[contig_full_grouping$clonotype_IGH == "none, dual IGH" | contig_full_grouping$clonotype_IGH == "none, only one chain" ] <- "none"
  contig_full_grouping$clonotype_IGKL[contig_full_grouping$clonotype_IGKL == "none, dual IGH" | contig_full_grouping$clonotype_IGKL == "none, only one chain" ] <- "none"

  #First, remove NA values from the clonotype numbers
  test <- contig_full_grouping %>%
    unite(clonotype_prelim, clonotype_IGH, clonotype_IGKL, na.rm = T, sep = ",")
  test$clonotype_prelim[test$clonotype_prelim == "none,none"] <- "none"

  #Then, make a concatenated preliminary clonotype based on the IGH and IGKL clonotype
  test2 <- paste.data.frame(test, term = "clonotype_prelim", group = "cell_barcode", collapse = "_")

  contig <- left_join(contig, test2, by = "cell_barcode")

  #Since we are going to cluster based on IGH sequences, only keep IGH sequences
  cluster_contig <- contig[contig$chain == "IGH" & contig$clonotype_prelim != "none",]
  cluster_contig$clonotype_new <- NA

  #Split the cells according to their donor
  group <- split(cluster_contig, list(cluster_contig$donor), drop = TRUE)

  #Split each donor group according to their preliminary clonotypes
  group <- lapply(X = group, function(x) split(x, x$clonotype_prelim, drop = T))


  #Establish the clonotype numbers
  clonotypes <- c()
  max_clono <- nrow(cluster_contig) #the maximum amount of clonotypes that can be generated is the number of cells in the object

  normal_nt <- c("A", "T", "G", "C") #establish what normal nucleotides are

  free_clonotypes <- setdiff(seq(1:max_clono), clonotypes) #set the clonotype numbers that are not in use

  print("Running the clonotyping process, hang in there:")
  progress_bar = txtProgressBar(min=0, max=length(group), style = 1, char="=") #set up the progress bar

  for (g in 1:length(group)){ #run for each donor
    subgroup <- group[[g]]

    group_start <- contig[contig$donor == names(group)[[g]],] #save the starting contig
    groups_to_be_checked <- c()


    for (s in 1:length(subgroup)){
      free_clonotypes <- setdiff(seq(1:max_clono), clonotypes)


      if(nrow(subgroup[[s]]) < 2){
        subgroup[[s]]$clonotype_new <- free_clonotypes[[1]] #if there is only one cell in the IGH/IGKL group, then no clonotyping needs to be done: this cell is unique
        clonotypes <- c(clonotypes, subgroup[[s]]$clonotype_new) #add the clonotype number to the ones in use

      } else{
        groups_to_be_checked <- c(groups_to_be_checked, s) #save the groups larger than one cell: these need to be further clonotyped

      }
    }

    for (s in groups_to_be_checked){ #check every group with more than one entry

      free_clonotypes <- setdiff(seq(1:max_clono), clonotypes)

      investigation_group <- subgroup[[s]] #get the group out


      for (x in 1:nrow(investigation_group)){ #for each entry in the data frame

        if(length(unique(investigation_group$clonotype_new)) == 1 && !any(is.na(investigation_group$clonotype_new))){
          next #if all cells already have the same clonotype, it does not make sense to continue

        }

        string_split <- unlist(strsplit(investigation_group$cdr3_nt[[x]], "", fixed = TRUE)) #make a vector of the individual nucleotides

        for (z in 1:nrow(investigation_group)){ #compare against each other entry in the data frame

          if (x == z){
            next #do not compare a sequence against itself
          }

          if(length(unique(investigation_group$clonotype_new)) == 1 && !any(is.na(investigation_group$clonotype_new))){
            next #if all cells already have the same clonotype, it does not make sense to continue

          }

          other_split <- unlist(strsplit(investigation_group$cdr3_nt[[z]], "", fixed = TRUE)) #make a vector of the individual nucleotides

          min_length <- min(length(string_split), length(other_split)) #the length to be compared over is the length of the smallest string

          score <- 0 #starting score is zero

          if (levenshtein == F){
            for (n in seq(1:min_length)){ #for each nucleotide, up to the length of the shortest sting
              if (other_split[n] != string_split[n] && other_split[n] %fin% normal_nt && string_split[n] %fin% normal_nt){ #add 1 to the mutation score if the nucleotides do not match up
                score = score + 1
              } else{
                score = score #do not add to the score if the nucleotides do match up or a non-standard nucleotide was used
              }
            }

          } else if (levenshtein == T){
            score <- adist(investigation_group$cdr3_nt[[x]], investigation_group$cdr3_nt[[z]], ignore.case = F)
          }

          if (score > threshold){ #if the mutation score exceeds a threshold (more than the threshold number of nucleotide mismatches)
            next #leave the clonotype as is

          } else{ #if the mutation score is below the threshold

            if(is.na(investigation_group$clonotype_new[[x]]) && is.na(investigation_group$clonotype_new[[z]])){
              investigation_group$clonotype_new[[x]] <- free_clonotypes[[1]] #if they do match, but the input sequence has no clonotype number, assign the input cell the first unused clonotype number
              investigation_group$clonotype_new[[z]] <- investigation_group$clonotype_new[[x]] #assign the comparison cell the clonotype number of the input sequence

            } else if (!is.na(investigation_group$clonotype_new[[x]]) && is.na(investigation_group$clonotype_new[[z]])){
              investigation_group$clonotype_new[[z]] <- investigation_group$clonotype_new[[x]] #if they do match, and the input cell already has a clonotype number, assign the comparison cell the clonotype number of the input sequence

            } else if (is.na(investigation_group$clonotype_new[[x]]) && !is.na(investigation_group$clonotype_new[[z]])){
              investigation_group$clonotype_new[[x]] <- investigation_group$clonotype_new[[z]] #if they do match, and only the comparison cell has a clonotype number, assign the input cell the clonotype number of the comparison cell

            } else if (!is.na(investigation_group$clonotype_new[[x]]) && !is.na(investigation_group$clonotype_new[[z]]) && investigation_group$clonotype_new[[x]] == investigation_group$clonotype_new[[z]]){ #if both cells already have a clonotype number
              next #if both sequences already have matching clonotype numbers, do not change anything

            } else if (!is.na(investigation_group$clonotype_new[[x]]) && !is.na(investigation_group$clonotype_new[[z]]) && investigation_group$clonotype_new[[x]] != investigation_group$clonotype_new[[z]]){
              investigation_group$clonotype_new[investigation_group$clonotype_new == investigation_group$clonotype_new[[z]]] <- investigation_group$clonotype_new[[x]]
              investigation_group$clonotype_new[[z]] <- investigation_group$clonotype_new[[x]] #if both cells have different clonotype numbers, assign all cells with the clonotype of the comparison cell to the clonotype of the input cell (merge the clonotypes)
            }

            clonotypes <- c(clonotypes, investigation_group$clonotype_new) #update which clonotype numbers are in use
            free_clonotypes <- setdiff(seq(1:max_clono), clonotypes)

          } #end of clonotyping



        } #end of z loop


      } #end of x loop

      #Fill in clonotype numbers for all cells that did not get a clonotype number at the end of the clonotyping process
      investigation_group$clonotype_new[is.na(investigation_group$clonotype_new)] <- free_clonotypes[1:length(investigation_group$clonotype_new[is.na(investigation_group$clonotype_new)])] #give all remaining cells their own clonotype
      clonotypes <- c(clonotypes, investigation_group$clonotype_new)
      free_clonotypes <- setdiff(seq(1:max_clono), clonotypes)

      investigation_group <- as.data.table(investigation_group) #convert to data table

      subgroup[[s]] <- investigation_group #update the list with the clonotyped object


    } #end of s loop
    group[[g]]  <- rbindlist(subgroup) #get the clonotyped object per donor in a data table

    group[[g]]$clonotype_new <- as.character(group[[g]]$clonotype_new)

    #double check if the clonotyping went correctly
    check <- group[[g]] %>%
      group_by(clonotype_new, v_gene, d_gene, j_gene) %>%
      count()
    check$correct <- "yes"
    check$correct[duplicated(check$clonotype_new)] <- "no" #if a clonotype has two entries, it means that they use different genes and the assignment did not go correctly

    check_wrong <- check[check$correct == "no",]

    wrong_cells <- group[[g]][group[[g]]$clonotype_new %fin% check_wrong$clonotype_new,] #get all cells with the incorrect clonotypes out

    #cells with the same end clonotype, but different prelim clonotypes, need to be split:
    if (nrow(wrong_cells) > 0){ #if there are any cells wrongly classified
      for (i in 1:nrow(wrong_cells)){
        free_clonotypes <- setdiff(seq(1:max_clono), clonotypes)
        same_clonotype_cells <- wrong_cells[wrong_cells$clonotype_new == wrong_cells$clonotype_new[[i]],] #find cells with the same new clonotype
        if (length(unique(same_clonotype_cells$clonotype_prelim)) > 1){ #if the cells of the same new clonotype have different prelim clonotypes
          wrong_cells$clonotype_new[[i]] <- free_clonotypes[[1]] #assign the current row a new clonotype
          wrong_cells$clonotype_new[wrong_cells$clonotype_prelim == wrong_cells$clonotype_prelim[[i]]] <- wrong_cells$clonotype_new[[i]] #also change the clonotype for the cells that do match perfectly with the prelim of the current row
          clonotypes <- c(clonotypes, wrong_cells$clonotype_new)
        }
      }
    }



    #Add any missing heavy chain chain data back in
    group_missed <- group_start[!(group_start$contig_id %fin% group[[g]]$contig_id),]
    group_missed$clonotype_new <- NA


    if (nrow(group_missed) > 0){
      group[[g]] <- rbind(group[[g]], group_missed)
      group[[g]]$clonotype_new[group[[g]]$clonotype_prelim == "none"] <- "none"
    }

    #fill in the clonotype numbers of the light chain data based on the heavy chain clonotypes
    group[[g]] <-  group[[g]] %>%
      group_by(cell_barcode) %>%
      fill(clonotype_new) %>%
      fill(clonotype_new, .direction = "up") %>%
      ungroup() %>%
      as.data.table()

    setTxtProgressBar(progress_bar, value = g) #add to the progress bar


  } #end of g loop

  print("Done clonotyping, wrapping up...")

  #Add the data of all donors back together
  group_data <- rbindlist(group)

  #Add more useful information to the "none" clonotype
  group_data$clonotype_new[group_data$clonotype_IGH == "none, dual IGH"] <- "none, dual IGH"
  group_data$clonotype_new[group_data$clonotype_IGH == "none, only one chain"] <- "none, only one chain"


  #Check if the output has the same amount of entries as the input:
  contig_missed <- contig_start[!(contig_start$contig_id %fin% group_data$contig_id),]

  if (nrow(contig_missed) > 0){ #if some entries were missed, add those back in
    contig_missed$clonotype_new <- NA
    contig_missed$clonotype_IGH <- "none, missing from end object"
    contig_missed$clonotype_IGKL <- "none, missing from end object"
    contig_missed$clonotype_prelim <- "none, missing from end object"

    group_data <- rbind(group_data, contig_missed)
  }

  #fill in the clonotype numbers of the missing light chain data based on the heavy chain clonotypes
  group_data <-  group_data %>%
    group_by(cell_barcode) %>%
    fill(clonotype_new) %>%
    fill(clonotype_new, .direction = "up") %>%
    ungroup() %>%
    as.data.table()

  group_data$clonotype_new[group_data$cdr3_nt == ""] <- "none, no CDR3" #assign clonotypes to all entries with no CDR3 data


  #Just in case there are NA values left:
  group_data$clonotype_new[is.na(group_data$clonotype_new)] <- "none, were not assigned clonotype"


  group_data$clonotype_IGH <- NULL
  group_data$clonotype_IGKL <- NULL
  group_data$clonotype_prelim <- NULL

  group_data <- as.data.frame(group_data)

  #Close the progress bar
  print("Finished! Returning the object to you now..")
  close(progress_bar)

  return(group_data)
}


#' Test whether the clonotype() function ran correctly
#'
#' This function tests whether the clonotypes generated by clonotype() really have identical heavy chain VDJ assignments and a maximum of nucleotide mismatches in the CDR3 given by the threshold.
#' It also shows how many cells were not assigned a clonotype, and how many cellss have multiple light chain assignments.
#' @param contig The clonotyped VDJ contigs, which is the result of the clonotype() function..
#' @param threshold The allowed number of mismatches (as an integer) in the CDR3 region between two sequences in a clonotype. This should be the same threshold as used in the clonotype() function.
#' @param levenshtein Should the levenshtein distance be used in the checks? Set this to be equal to the clonotype() function. Default is FALSE.
#' @keywords checks
#' @returns Prints out the amount of non-clonotyped cells, the number of cells with multiple light chain assignments, and whether the object passed the heavy chain and CDR3 checks.
#' @export
#' @examples contig <- clonotype(contig, threshold = 2, levenshtein = F)
#'            contig <- test.result(contig, threshold = 2, levenshtein = F)

test.result <- function(contig, threshold, levenshtein = F){

  #Set back all "none"clonotype entries
  contig_none <- contig[grepl("none", contig$clonotype_new),]
  print("Sequences with no clonotype assignment: ")
  print(summary(as.factor(contig_none$clonotype_new)))

  contig$clonotype_new[grepl("none", contig$clonotype_new)] <- "none"


  #IGH chain check
  check_IGH <- contig %>%
    filter(chain == "IGH") %>%
    filter(clonotype_new != "none") %>%
    group_by(clonotype_new, v_gene, d_gene, j_gene) %>%
    dplyr::count()

  dup1 <- check_IGH[duplicated(check_IGH$clonotype_new),]
  if (nrow(dup1) > 0){
    passed_IGH <- F
  } else{
    passed_IGH <- T
  }

  check_IGKL <- contig %>%
    filter(chain == "IGK" | chain == "IGL") %>%
    filter(clonotype_new != "none") %>%
    group_by(clonotype_new, v_gene, j_gene) %>%
    dplyr::count()

  dup2 <- check_IGKL[duplicated(check_IGKL$clonotype_new),]
  #does not have to be zero due to double IGKL
  print(paste0(nrow(dup2), " clonotypes with more than one possible IGKL chain were found."))

  check_CDR3 <- contig %>%
    filter(chain == "IGH") %>%
    filter(clonotype_new != "none") %>%
    group_by(clonotype_new, v_gene, d_gene, j_gene, cdr3_nt) %>%
    dplyr::count()

  differing_cdr3 <- check_CDR3[duplicated(check_CDR3$clonotype_new),]
  differing_cdr3 <- check_CDR3[check_CDR3$clonotype_new %fin% differing_cdr3$clonotype_new,]

  check_distances <- function(differing_cdr3, threshold, levenshtein){
    results_data <- data.frame(matrix(ncol = 1, nrow = length(unique(differing_cdr3$clonotype_new))))
    colnames(results_data) <- "cleared?"
    rownames(results_data) <- unique(differing_cdr3$clonotype_new)

    for (i in unique(differing_cdr3$clonotype_new)){
      group <- differing_cdr3[differing_cdr3$clonotype_new == i,]
      if (nrow(group) == 2){
        if (levenshtein == T){
          distance <- adist(group$cdr3_nt[[1]],  group$cdr3_nt[[2]], ignore.case = F)
        } else{
          distance <- adist(group$cdr3_nt[[1]], group$cdr3_nt[[2]], costs = c(ins = 100, del = 100, sub = 1)) #costs of indels are set to 100 to prevent these from happening

        }
        if (distance > threshold){
          results_data[i,] <- "no"
        } else{
          results_data[i,] <- "yes"

        }
      }

      if (nrow(group) > 2){
        distances <- c()
        group_result <- c()
        for (x in 1:nrow(group)){
          for (z in 1:nrow(group)){
            if (x==z){
              next
            }
            if (levenshtein == T){
              distance <- adist(group$cdr3_nt[[x]],  group$cdr3_nt[[z]], ignore.case = F)
            } else{
              distance <- adist(group$cdr3_nt[[x]], group$cdr3_nt[[z]], costs = c(ins = 100, del = 100, sub = 1)) #costs of indels are set to 100 to prevent these from happening

            }
            distances <- c(distances, distance)

          }


          if (any(distances <= threshold)){ #if a group member has a distance below the threshold to any other rmember, it passes checks
            group_result[[i]] <- "yes"
          } else{
            group_result[[i]] <- "no"

          }

        }

        if (any(group_result == "no")){
          results_data[i,] <- "no"

        } else{
          results_data[i,] <- "yes"

        }

      }
    }
    return(results_data)
  }

  check_CDR3.res <- check_distances(differing_cdr3, threshold = threshold, levenshtein = levenshtein)
  failed.checks <- rownames(check_CDR3.res)[which(check_CDR3.res == "no")] #these need to be zero

  failed.checks <- check_CDR3[check_CDR3$clonotype_new %fin% failed.checks,]
  if (nrow(failed.checks) > 0){
    passed_CDR3 <- F
  } else{
    passed_CDR3 <- T
  }


  if (passed_IGH == F && passed_CDR3 == F){
    print("Result: object did not meet IGH and CDR3 criteria.")
  } else if (passed_IGH == F && passed_CDR3 == T){
    print("Result: object did not meet IGH criteria.")

  } else if (passed_IGH == T && passed_CDR3 == F){
    print("Result: object did not meet CDR3 criteria. Did you check if you entered the correct threshold?")

  } else if (passed_IGH == T && passed_CDR3 == T){
    print("Result: object passed all tests.")

  }

  return(contig)

}


#' Assign autoreactivity
#'
#' This function assigns cells as being autoreactive or not based on the CDR3 amino acid sequences of known autoreactive cells.
#' If the CDR3 aa sequence of a cell matches a known autoreactive sequence, the cell and its entire clonotype are annotated as autoreactive.
#' @param contig The clonotyped and tested VDJ contigs as is the result from clonotype() and test.result().
#' @param sequences The amino acid sequences of B cells that are known to be autoreactive. Format is data.frame. This data.frame should also contain which donor the sequence was found in, and a name for the specific reactivity (for example, CARAAA-Ro52)
#' @keywords autoreactivity
#' @returns The original contig file, with two new columns: their reactivities (named the same as in the sequences data.frame) and their antigenicity (for example, RF or La).
#' @export
#' @examples contig <- clonotype(contig, threshold = 2, levenshtein = F)
#'             contig <- assign_autoreactivity(contig, sequences)

assign_autoreactivity <- function(contig, sequences){

  if(all(c("reactivity", "antigen") %fin% colnames(contig))){
    print("Warning: overriding previous reactivity and antigen assignments")
  }

  contig$reactivity <- NA
  contig$antigen <- NA
  antigen_vector <- c()


  if(!(all(c("cdr3", "donor", "reactivity") %in% colnames(sequences)))){
    stop("Please add columns in the sequences input stating the cdr3, donor and reactivity, and name them accordingly.")
  }



  for (i in 1:nrow(sequences)){
    for (j in 1:nrow(contig)){
      if (sequences$cdr3[[i]] == contig$cdr3[[j]] && (sequences$donor[[i]] == contig$donor[[j]] | contig$donor[[j]] == "AIM")
          && contig$clonotype_new[[j]] != "none" && contig$clonotype_new[[j]] != "none, dual IGH" && contig$clonotype_new[[j]] != "none, only one chain"){ #if the donor and autoreactive sequence match, or the donor is AIM
        if (is.na(contig$reactivity[[j]])){
          contig$reactivity[[j]] <- sequences$reactivity[[i]] #mark the cell as autoreactive if it matches the sequence
        } else if (!grepl(sequences$reactivity[[i]], contig$reactivity[[j]])){
          contig$reactivity[[j]] <- paste0(contig$reactivity[[j]], "," ,sequences$reactivity[[i]]) #if the cell already had a reactivity, but not the one being tested against, add the new reactivity
        } else{
          next #if the cell already has that reactivity, skip
        }
        contig$reactivity[contig$clonotype_new == contig$clonotype[[j]]] <- contig$reactivity[[j]] #assign cells with the same clonotype the same reactivity
      } else{
        next
      }
    }
  }

  #now, assign the antigenicity:
  for (i in 1:nrow(contig)){
    row <- contig[i,]
    antigen <- NA
    if (grepl("-RF", row$reactivity, fixed = TRUE)){
      antigen <- "RF"
    }

    if (grepl("-Ro60", row$reactivity, fixed = TRUE)){
      if (is.na(antigen)){
        antigen <- "Ro60"

      } else{
        antigen <- paste0(antigen, ",", "Ro60")
      }
    }

    if (grepl("-Ro52", row$reactivity, fixed = TRUE)){
      if (is.na(antigen)){
        antigen <- "Ro52"

      } else{
        antigen <- paste0(antigen, ",", "Ro52")
      }
    }

    if (grepl("-La", row$reactivity, fixed = TRUE)){
      if (is.na(antigen)){
        antigen <- "La"

      } else{
        antigen <- paste0(antigen, ",", "La")
      }
    }

    antigen_vector[i] <- antigen
  }

  contig$antigen <- antigen_vector


  return(contig)
}

#' Remove VDJ doublets
#'
#' This function removes doublets from gene expression matrices (and CITE-seq matrices, if applicable) based on the VDJ data of that matrix.
#' Dual IGH or TRB chains are very rare in lymphocytes and are likely the result of a doublet. Also, cells should either express a BCR or a TCR, not both.
#' This function removes cells with dual IGH or TRB chains, or cells that appear in both the VDJ B and VDJ T data from the gene expression and CITE count matrices.

#' @param matrix A matrix containing gene expression counts, or a list consisting of a gene expression count matrix and a CITE-seq count matrix. This function is directly applicable to 10X expression matrices generated by cellranger.
#' @param vdj_b A data.frame containing the VDJ B contig data
#' @param vdj_t A data.frame containing the VDJ T contig data
#' @keywords doublet, vdj
#' @returns The original count matrices with the doublet cells removed.
#' @export
#' @examples matrix <- remove_vdj_doublets(matrix, vdj_b, vdj_t)

remove_vdj_doublets <- function(matrix, vdj_b, vdj_t){
  stopifnot(require(dplyr))
  stopifnot(require(Matrix))

  doublets <- c()
  IGH_doublets <- vdj_b %>%
    group_by(barcode) %>%
    filter(chain == "IGH") %>%
    dplyr::count()
  IGH_doublets <- IGH_doublets[IGH_doublets$n > 1,]
  IGH_doublets <- IGH_doublets$barcode
  print(paste0(length(IGH_doublets), " cells found with more than one IGH chain."))

  TRB_doublets <- vdj_t %>%
    group_by(barcode) %>%
    filter(chain == "TRB") %>%
    dplyr::count()
  TRB_doublets <- TRB_doublets[TRB_doublets$n > 1,]
  TRB_doublets <- TRB_doublets$barcode
  print(paste0(length(TRB_doublets), " cells found with more than one TRB chain."))

  T_B_doublets <- intersect(vdj_b$barcode, vdj_t$barcode) #find cells that have both a TCR and BCR entry
  doublets <- c(IGH_doublets, TRB_doublets, T_B_doublets)
  doublets <- doublets[!duplicated(doublets)]
  print(paste0(length(T_B_doublets), " cells found with a TCR and BCR."))


  #filter the gene expression (and CITE, if applicable)
  if (type(matrix) == "list"){
    matrix[["Gene Expression"]] <- matrix[["Gene Expression"]][,!(colnames(matrix[["Gene Expression"]]) %in% doublets)]
    matrix[["Antibody Capture"]] <- matrix[["Antibody Capture"]][,!(colnames(matrix[["Antibody Capture"]]) %in% doublets)]

  } else{
    matrix <- matrix[,!(colnames(matrix) %in% doublets)]

  }
  print(paste0(length(doublets), " doublets removed."))

  return(matrix)
}

#' Generate seurat objects from a list
#'
#' This is a simple function to generate a list of seurat objects from a list of count matrices containing gene expression data or gene expression data and CITE-seq data.
#' @param expression_matrices A named list of expression matrices. These expression matrices can be, for example, cellranger outputs.
#' @param min.cells The minimum amount of cells a gene needs to be expressed in to be kept in the seurat object
#' @param min.genes The minimum amount of genes a cell needs to express to be kept in the seurat object
#' @keywords seurat
#' @returns A named list of Seurat objects with RNA and CITE assays (wherever applicable)
#' @export
#' @examples list <- list(expression_matrix1, expression_matrix2)
#'           names(list) <- c("object1", "object2")
#'           seurat_list <- generate_seurat_objects(list)
#'
generate_seurat_objects <- function(expression_matrices, min.cells = 3, min.genes = 200){
  stopifnot(require(Seurat))

  seurat_objects <- list()
  for (i in seq_along(expression_matrices)) {
    #Create a new seurat object
    if (type(expression_matrices[[i]]) == "list"){
      invisible(seurat_object <- CreateSeuratObject(counts = expression_matrices[[i]][["Gene Expression"]], min.cells = min.cells, min.genes = min.genes, project = names(expression_matrices)[[i]]))
    } else{
      invisible(seurat_object <- CreateSeuratObject(counts = expression_matrices[[i]], min.cells = min.cells, min.genes = min.genes, project = names(expression_matrices)[[i]]))

    }
    seurat_object@meta.data$barcode <- rownames(seurat_object@meta.data)

    #Set the condition and donor code
    seurat_object@meta.data$condition <- "pSS"
    seurat_object@meta.data$donor <- names(expression_matrices)[[i]]

    if (type(expression_matrices[[i]]) == "list"){
      if ("Antibody Capture" %in% names(expression_matrices[[i]])){
        cite <- expression_matrices[[i]][["Antibody Capture"]]
        suppressWarnings(seurat_object[["CITE"]] <- CreateAssayObject(cite))


      } else {
        suppressWarnings(rm(cite))
      }

    }
    #add to the list
    seurat_objects[[i]] <- seurat_object
    names(seurat_objects)[[i]] <- names(expression_matrices)[[i]]


  }
  return(seurat_objects)
}


#' Add VDJ data to seurat objects in a list
#'
#' This function adds VDJ data to a list of seurat objects using scRepertoire. The clonecall is set to "strict".
#' @param list A named list of seurat objects. This can be generated by the generate_seurat_objects() function.
#' @param contigs A named list of contigs from VDJ data. This list should be in the same order as the list of seurat objects.
#' @keywords seurat, screpertoire, vdj
#' @returns A named list of Seurat objects with VDJ data added.
#' @export
#' @examples list <- list(expression_matrix1, expression_matrix2)
#'           names(list) <- c("object1", "object2")
#'           seurat_list <- generate_seurat_objects(list)
#'          contig_list <- list(contig_1, contig_2)
#'          names(contig_list) <- c("object1", "object2")
#'           seurat_list <- add_BCR(seurat_list, contig_list)

add_BCR <- function(list, contigs){
  stopifnot(require(scRepertoire))
  stopifnot(require(Seurat))

  for (i in seq_along(list)){
    if (length(list) == length(contigs)){
      BCR <- scRepertoire::combineBCR(contigs[[i]], samples = c("1"), ID =c("1"), removeNA = TRUE, removeMulti = FALSE)
      BCR <- BCR[[1]]
      BCR$barcode <- gsub("1_1_", "", BCR$barcode)
      seurat_object <- list[[i]]

      seurat_object <- combineExpression.new(df = BCR,sc =  seurat_object,
                                             cloneCall="strict",
                                             proportion = FALSE,
                                             cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded = 500))

      seurat_object@meta.data$barcode <- rownames(seurat_object@meta.data)
      seurat_object@meta.data$cloneType <- as.factor(seurat_object@meta.data$cloneType)
      name <- names(list)[[i]]
      list[[i]] <- seurat_object
      names(list)[[i]] <- name
    } else {
      print("Please supply as many contig information tables as seurat objects")
    }
  }
  return(list)
}


#' Remove contaminating genes
#'
#' This function removes contaminating genes and genes that are not of interest (mitochondrial, ribosomal, MALAT1, GAPDH, microRNAs, ncRNAs) from the RNA assays of a list of seurat objects.
#' This function should be run prior to normalisation.
#' @param list A named list of seurat objects. This can be generated by the generate_seurat_objects() function.
#' @param seurat_list A named list of seurat objects, which can be generated by the generate_seurat_objects() function.
#' @keywords seurat, RNA
#' @returns A named list of Seurat objects with contaminating and non-interesting genes removed.
#' @export
#' @examples list <- list(expression_matrix1, expression_matrix2)
#'           names(list) <- c("object1", "object2")
#'           seurat_list <- generate_seurat_objects(list)
#'           seurat_list <- remove_genes(seurat_list)

remove_genes <- function(seurat_list){
  stopifnot(require(Seurat))
  seurat_list_filtered <- list()
  suppressWarnings(rm(cite))

  for (i in seq_along(seurat_list)){
    suppressWarnings(rm(cite))
    if ("CITE" %in% names(seurat_list[[i]])){
      cite <- seurat_list[[i]][["CITE"]]
    } else {
      suppressWarnings(rm(cite))
    }

    rem_genes <- grep("^mt-|^RPS|^MRPS|^RPL|^MRPL|MALAT1|^MIR|GAPDH|^RP|\\.+", rownames(seurat_list[[i]][["RNA"]]@counts) , ignore.case=T, value=T)
    #finds these genes:
    #MT- (Mitochondrial proteins),
    #RP- RPS- MRPS- RPL- MRPL- ((mitochondrial)Ribosomal proteins)
    #Contaminating genes like MALAT1 and GAPDH
    #non-coding RNAs (these have a dot in their gene name)

    counts <- seurat_list[[i]][["RNA"]]@counts #extracts the counts table
    counts <- counts[-(which(rownames(counts) %in% rem_genes)),] #looks for which genes are not in the list of genes to be removed
    seurat_list_filtered[[i]] <- subset(seurat_list[[i]], features = rownames(counts)) #subsets the seurat based on the genes we want to keep
    names(seurat_list_filtered)[[i]] <- names(seurat_list)[[i]]
    rm(counts)

    if (exists("cite", mode = "S4")){
      seurat_list_filtered[[i]][["CITE"]] <- cite
    }
  }
  return(seurat_list_filtered)
}

#' Annotate GO and KEGG pathways of all clusters
#'
#' This function performs GO and KEGG pathway analysis (using clusterProfiler) for each cluster based on the marker genes found with seurats FindAllMarkers() function. It does not discriminate between upregulated and downregulated markers.
#' If this distinction is desired, the data frame should be split and the function should be ran on each dataframe.
#' If gene names appear as their non-standard names, the function will try to convert these to their standard names.

#' @param list A named list of seurat objects. This can be generated by the generate_seurat_objects() function.
#' @param seurat_list A named list of seurat objects, which can be generated by the generate_seurat_objects() function.
#' @keywords GO, KEGG, markers
#' @returns A list with a GO and a KEGG result for each cluster. The dotplot() function of clusterProfiler can be used on each list element.
#' @export
#' @examples markers <- FindAllMarkers(seurat)
#'           cluster_annotation <- annotate_clusters(markers)
#'           dotplot(cluster_annotation[["GO"]])

annotate_clusters <- function(markers){
  stopifnot(require(clusterProfiler))
  stopifnot(require(limma))

  gene_entrez <-  NULL
  gene_symbol <- NULL

  if (!("annotation" %in% colnames(markers))){
    markers$annotation <- markers$cluster
  }

  levels <- unique(markers$annotation)

  for (i in levels){
    genes <- subset(markers, markers$annotation == i) #retrieve the marker data for a specific cluster
    genes <- unlist(genes$gene)
    tryCatch({
      suppressMessages(genes_convert <- bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db"))

    }, error = function(e){genes_convert <<- c()})

    if (length(genes_convert) > 0){
      genes_unconvert <- genes[!(genes %in% genes_convert$SYMBOL)]
      genes_convert <- genes_convert$ENTREZID
    } else {
      genes_unconvert <- genes
      genes_convert <- NA
    }


    if (length(genes_unconvert) > 0){

      genes_unconvert_map <- limma::alias2Symbol(genes_unconvert, species = "Hs", expand.symbols = FALSE) #some genes do not convert because the name is an alias for the official name. Map them to their official name.
      #attempt to convert the genes to ENTREZID if their symbol was converted:
      tryCatch({
        suppressMessages(genes_unconvert_map <- bitr(genes_unconvert_map, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db"))

      }, error = function(e){print(paste0("Error: these genes could not be converted on second try: ", genes_unconvert))})
      genes_unconvert.2 <- genes_unconvert[!(genes_unconvert %in% genes_unconvert_map$SYMBOL)]

      print(paste0("Not converted on second try: ", genes_unconvert.2))

      genes_unconvert_map <- genes_unconvert_map$ENTREZID

      genes_convert <- c(genes_convert, genes_unconvert_map)
    }


    gene_entrez[[length(gene_entrez)+1]] <- genes_convert #add the genes of that cluster to the list
    names(gene_entrez)[[length(gene_entrez)]] <-i
    gene_symbol[[length(gene_entrez) +1]] <- genes
    names(gene_symbol)[[length(gene_entrez)]] <- i


  }

  ck <- compareCluster(geneCluster = gene_entrez, fun = enrichGO, OrgDb = org.Hs.eg.db, pvalueCutoff = 0.05, qvalueCutoff = 0.05) #perform the comparison
  summary_ck <- as.data.frame(ck)


  #compare KEGG pathways between clusters
  cke <- compareCluster(geneCluster = gene_entrez, fun = enrichKEGG, pvalueCutoff = 0.05, qvalueCutoff = 0.05) #perform the comparison
  summary_cke <- as.data.frame(cke)

  save <- list(ck, cke)
  names(save) <- c("GO", "KEGG")

  return(save)

}

#' Group Azimuth l2 annotations
#'
#' This function groups Azimuth level 2 annotations into larger, overarching groups. This is mainly meant for seurat objects consisting of of B cells.
#' This function should be run prior to normalisation.
#' @param seurat_object A seurat object on which RunAzimuth() has been called.
#' @keywords seurat, Azimuth
#' @returns A new column in the metadata of the seurat object, called top_level_annot, which contains the grouped annotations.
#' @export
#' @examples seurat <- RunAzimuth(seurat, ref = "tonsilref")
#'           seurat <- add_top_level(seurat)
#'
add_top_level <- function(seurat_object){

  stopifnot(require(Seurat))
  stopifnot(require(dplyr))

  if (!("predicted.celltype.l2" %in% colnames(seurat_object@meta.data))){
    stop("Run Azimuth on the seurat object first.")
  }

  annotations_azimuth <- seurat_object@meta.data %>%
    group_by(predicted.celltype.l1, predicted.celltype.l2) %>%
    dplyr::count()

  annotations_azimuth$top_level_annot <- NA

  for (i in 1:nrow(annotations_azimuth)){

    annotation <- annotations_azimuth$predicted.celltype.l2[[i]]
    if (grepl("Mast|macrophage|monocyte|dendritic|DC", annotation, ignore.case = F)){
      annotations_azimuth$top_level_annot[[i]] <- "Myeloid"
    } else if (grepl("NK|ILC", annotation, ignore.case = F)){
      annotations_azimuth$top_level_annot[[i]] <- "NK/ILC"
    } else if (grepl("NK|ILC", annotation, ignore.case = F)){
      annotations_azimuth$top_level_annot[[i]] <- "NK/ILC"
    } else if (grepl("T|Naive", annotation, ignore.case = F)){
      annotations_azimuth$top_level_annot[[i]] <- "T cell"
    } else if (grepl("Crypt", annotation, ignore.case = F)){
      annotations_azimuth$top_level_annot[[i]] <- "Epithelial"


    } else if (grepl("NBC", annotation, ignore.case = F)){
      annotations_azimuth$top_level_annot[[i]] <- "NBC" #annotate naive B cells

    } else if (grepl("MBC", annotation, ignore.case = F) && !grepl("plasma|Plasma|PC", annotation, ignore.case = F)){
      annotations_azimuth$top_level_annot[[i]] <- "MBC" #annotate memory B cells

    } else if (grepl("MBC", annotation, ignore.case = F) && grepl("plasma|Plasma|PC", annotation, ignore.case = F)){
      annotations_azimuth$top_level_annot[[i]] <- "MBC-derived PC" #annotate memory B cell-derived plasma cells

    } else if (grepl("plasma|Plasma|PC", annotation, ignore.case = F) && grepl("IgG|IgA", annotation, ignore.case = F)){
      if (grepl("precursor|preMature|early", annotation, ignore.case = T)){
        annotations_azimuth$top_level_annot[[i]] <- "IgG/IgA precursor" #annotate mature plasma cells as precursor

      } else {
        annotations_azimuth$top_level_annot[[i]] <- "IgG/IgA" #annotate mature plasma cells as mature

      }

    } else if (grepl("plasma|Plasma|PC", annotation, ignore.case = F) && grepl("IgM|IgD", annotation, ignore.case = F)){
      if (grepl("precursor|preMature|early", annotation, ignore.case = T)){
        annotations_azimuth$top_level_annot[[i]] <- "IgM/IgD precursor" #annotate immature plasma cells as precursor

      } else{
        annotations_azimuth$top_level_annot[[i]] <- "IgM/IgD" ##annotate immature plasma cells as developed
      }

    } else if (grepl("DZ|Dark Zone", annotation, ignore.case = F) && grepl("LZ", annotation, ignore.case = F)){
      annotations_azimuth$top_level_annot[[i]] <- "DZ/LZ"

    } else if (grepl("DZ", annotation, ignore.case = F)){
      annotations_azimuth$top_level_annot[[i]] <- "centroblast"

    } else if (grepl("LZ", annotation, ignore.case = F)){
      annotations_azimuth$top_level_annot[[i]] <- "centrocyte"

    } else if (grepl("GC", annotation, ignore.case = F)){
      annotations_azimuth$top_level_annot[[i]] <- "GC precursor"

    } else {
      annotations_azimuth$top_level_annot[[i]] <- annotations_azimuth$predicted.celltype.l2[[i]]
    }



  }
  #take the top level annotations and ungroup the data frame
  annotations_azimuth <- annotations_azimuth[,c("predicted.celltype.l2", "top_level_annot")]
  annotations_azimuth <- as.data.frame(ungroup(annotations_azimuth))

  #Get the seurat metadata out
  seurat_annot <- seurat_object@meta.data

  #match the top level annotations to the level 2 annotations in the seurat metadata
  annotations_azimuth.full <- as.data.frame(annotations_azimuth$top_level_annot[match(seurat_annot$predicted.celltype.l2, annotations_azimuth$predicted.celltype.l2)])
  colnames(annotations_azimuth.full) <- c("top_level_annot")

  #add the top level annotations to the seurat metadata
  seurat_annot <- cbind(seurat_annot, annotations_azimuth.full)

  seurat_object@meta.data <- seurat_annot

  return(seurat_object)

}

#' Add dotplot annotation
#'
#' This function adds an annotation that can be used to separate groups in a dotplot. Cells with no autoreactivity assigned will be annotated based on their isotype,
#' and cells that are assigned an autoreactivity based on their group (RF, ANA or both). The amount of cells that are present within a group(cluster+ isotype/reactivity)
#' @param seurat A seurat object with the following columns in its metadata: annotation: the cluster names or grouped_annotation: the grouped annotation generated by group_annot(),
#' group: the group the cells belong to (RF, ANA, both or NA), c_gene: the isotype as predicted by cellranger or a similar program
#' @param mode Which annotation to use: the cluster annotation ("cluster") or the grouped annotation ("group") generated by the group_annot() function. Default is cluster.
#' @param include_IGHA_IGHD logical indicating whether cells with an IGHA or IGHD isotype should be annotated too.
#' @keywords seurat, dotplot
#' @returns A new column in the seurat metadata (dotplot_annot_res) giving the new annotation.
#' @export
#' @examples seurat <- add_dotplot_annot(seurat, mode = "cluster", include_IGHA_IGHD = T)
#'
add_dotplot_annot <- function(seurat, mode = "cluster", include_IGHA_IGHD = F){
  stopifnot(require(Seurat))
  stopifnot(require(dplyr))

  if(!("annotation" %in% colnames(seurat@meta.data)) & mode == "cluster"){
    stop("Annotation column is missing")
  }

  if(!("grouped_annotation" %in% colnames(seurat@meta.data)) & mode == "group"){
    stop("Grouped annotation column is missing")
  }

  if(!("c_gene" %in% colnames(seurat@meta.data))){
    stop("Isotype assignment is missing")
  }
  seurat@meta.data$dotplot_annot <- NA
  seurat@meta.data$n <- NULL

  if(mode == "cluster"){
    for (i in 1:nrow(seurat@meta.data)){
      if (!is.na(seurat@meta.data$group[[i]])){
        if (seurat@meta.data$group[[i]] %in% c("RF", "both")){
          seurat@meta.data$dotplot_annot[[i]] <- paste0(seurat@meta.data$annotation[[i]], "_RF")
        } else if (seurat@meta.data$group[[i]] == "ANA"){
          seurat@meta.data$dotplot_annot[[i]] <- paste0(seurat@meta.data$annotation[[i]], "_ANA")
        }

      } else{
        if (is.na(seurat@meta.data$c_gene[[i]])){
          next #cannot assign annotation without a c gene
        }
        if (seurat@meta.data$c_gene[[i]] == "IGHM"){
          seurat@meta.data$dotplot_annot[[i]] <- paste0(seurat@meta.data$annotation[[i]], "_IGHM")

        } else if (seurat@meta.data$c_gene[[i]] %in% c("IGHG1", "IGHG2", "IGHG3", "IGHG4")){
          seurat@meta.data$dotplot_annot[[i]] <- paste0(seurat@meta.data$annotation[[i]], "_IGHG")

        }
        if (include_IGHA_IGHD == T){
          if (seurat@meta.data$c_gene[[i]] == "IGHD"){
            seurat@meta.data$dotplot_annot[[i]] <- paste0(seurat@meta.data$annotation[[i]], "_IGHD")

          } else if (seurat@meta.data$c_gene[[i]] %in% c("IGHA1", "IGHA2")){
            seurat@meta.data$dotplot_annot[[i]] <- paste0(seurat@meta.data$annotation[[i]], "_IGHA")

          }
        }
      }

    }
  } else if (mode == "group"){
    for (i in 1:nrow(seurat@meta.data)){
      if (!is.na(seurat@meta.data$group[[i]])){
        if (seurat@meta.data$group[[i]] %in% c("RF", "both")){
          seurat@meta.data$dotplot_annot[[i]] <- paste0(seurat@meta.data$grouped_annotation[[i]], "_RF")
        } else if (seurat@meta.data$group[[i]] == "ANA"){
          seurat@meta.data$dotplot_annot[[i]] <- paste0(seurat@meta.data$grouped_annotation[[i]], "_ANA")
        }

      } else{
        if (is.na(seurat@meta.data$c_gene[[i]])){
          next #cannot assign annotation without a c gene
        }
        if (seurat@meta.data$c_gene[[i]] == "IGHM"){
          seurat@meta.data$dotplot_annot[[i]] <- paste0(seurat@meta.data$grouped_annotation[[i]], "_IGHM")

        } else if (seurat@meta.data$c_gene[[i]] %in% c("IGHG1", "IGHG2", "IGHG3", "IGHG4")){
          seurat@meta.data$dotplot_annot[[i]] <- paste0(seurat@meta.data$grouped_annotation[[i]], "_IGHG")

        }
        if (include_IGHA_IGHD == T){
          if (seurat@meta.data$c_gene[[i]] == "IGHD"){
            seurat@meta.data$dotplot_annot[[i]] <- paste0(seurat@meta.data$grouped_annotation[[i]], "_IGHD")

          } else if (seurat@meta.data$c_gene[[i]] %in% c("IGHA1", "IGHA2")){
            seurat@meta.data$dotplot_annot[[i]] <- paste0(seurat@meta.data$grouped_annotation[[i]], "_IGHA")

          }
        }
      }

    }
  } else {
    stop("Please select one of the following modes: cluster or group")
  }

  #add group sizes
  group_sizes <- seurat@meta.data %>%
    group_by(dotplot_annot) %>%
    filter(!is.na(dotplot_annot)) %>%
    dplyr::count()

  seurat@meta.data <- left_join(seurat@meta.data, group_sizes, by = "dotplot_annot")
  rownames(seurat@meta.data) <- seurat@meta.data$cell_barcode

  seurat@meta.data$dotplot_annot_res <- paste0(seurat@meta.data$dotplot_annot, " (", seurat@meta.data$n, ")")
  seurat@meta.data$dotplot_annot_res[is.na(seurat@meta.data$dotplot_annot)] <- "none"

  seurat@meta.data$dotplot_annot <- NULL
  seurat@meta.data$n <- NULL


  return(seurat)
}

#' Group clusters into larger groups
#'
#' Add an annotation that groups your Seurat clusters into 3 groups (MBC, PC, GC) per tissue (SG, LN, PBMC). Any cluster that does not belong to one of the three groups is assigned "Other".
#' The function throws an error if a cluster is assigned to multiple groups, so make sure no cluster number appears in multiple vectors. The cluster numbering starts counting from zero.
#' @param seurat A seurat object with the following columns in its metadata: annotation: the cluster names, tissue: the tissue of origin of that cell
#' @param annotation A vector with the unique annotation for each cluster. This can, for example, be the vector you use to rename the seurat clusters.
#' @param gc_index A vector of the cluster number of all clusters belonging to the germinal center (GC) class.
#' @param mbc_index A vector of the cluster number of all clusters belonging to the memory B cell (MBC) class.
#' @param plasma_index A vector of the cluster number of all clusters belonging to the plasma cell (PC) class.
#' @param nbc_index A vector of the cluster number of all clusters belonging to the naive B cell cell (NMC) class.
#' @keywords seurat, grouping, cluster
#' @returns A new column called grouped_annotation in the seurat metadata giving the new annotation.
#' @export
#' @examples annotation <- c("Germinal center", "CD27+ MBC", "CXCR4+ MBC", "IGHM+ PC", "IGHG+ PC")
#' seurat <- group_annot(seurat, annotation, gc_index = c(0), mbc_index = c(1,2), plasma_index = c(3,4))
#'
group_annot <- function(seurat, annotation, mbc_index, gc_index, plasma_index, nbc_index){
  stopifnot(require(Seurat))
  seurat@meta.data$grouped_annotation <- NA

  intersect_mbc_gc <- intersect(mbc_index, gc_index)
  intersect_mbc_plasma <- intersect(mbc_index, plasma_index)
  intersect_gc_plasma <- intersect(plasma_index, gc_index)
  intersect_mbc_naive <- intersect(mbc_index, nbc_index)
  intersect_plasma_naive <- intersect(plasma_index, nbc_index)
  intersect_gc_naive <- intersect(gc_index, nbc_index)


  if (length(intersect_gc_plasma) > 0 | length(intersect_mbc_gc) > 0 | length(intersect_mbc_plasma) > 0 |
      intersect_mbc_naive > 0 | intersect_plasma_naive > 0 | intersect_gc_naive > 0){
    error_intersect <- c(intersect_gc_plasma, intersect_mbc_gc, intersect_mbc_plasma, intersect_mbc_naive, intersect_plasma_naive, intersect_gc_naive)
    stop(paste0("Clusters cannot belong to two groups. Error occurred at cluster(s): ", unique(error_intersect)))
  }

  if (all(!(c("Lymph node", "Salivary gland", "Blood") %in% seurat@meta.data$tissue))){
    stop("None of the requested tissues were found in the seurat metadata: lymph node, salivary gland, blood")
  }

  mbc_index <- mbc_index + 1 #correct for the fact that a vector cannot start at 0
  gc_index <- gc_index + 1 #correct for the fact that a vector cannot start at 0
  plasma_index <- plasma_index + 1 #correct for the fact that a vector cannot start at 0
 nbc_index <- nbc_index + 1 #correct for the fact that a vector cannot start at 0

  #LN MBC, PC, GC (use one number higher than cluster number as annotation is a vector and cannot have element 0)
  seurat@meta.data$grouped_annotation[seurat@meta.data$annotation %in% annotation[mbc_index] & seurat@meta.data$tissue == "Lymph node"] <- "LN MBC"

  seurat@meta.data$grouped_annotation[seurat@meta.data$annotation %in% annotation[plasma_index] & seurat@meta.data$tissue == "Lymph node"] <- "LN Plasma"

  seurat@meta.data$grouped_annotation[seurat@meta.data$annotation %in% annotation[gc_index] & seurat@meta.data$tissue == "Lymph node"] <- "LN GC"

  seurat@meta.data$grouped_annotation[seurat@meta.data$annotation %in% annotation[nbc_index] & seurat@meta.data$tissue == "Lymph node"] <- "LN NBC"


  #SG MBC, PC, GC
  seurat@meta.data$grouped_annotation[seurat@meta.data$annotation %in%  annotation[mbc_index] & seurat@meta.data$tissue == "Salivary gland"] <- "SG MBC"

  seurat@meta.data$grouped_annotation[seurat@meta.data$annotation %in% annotation[plasma_index] & seurat@meta.data$tissue == "Salivary gland"] <- "SG Plasma"

  seurat@meta.data$grouped_annotation[seurat@meta.data$annotation %in% annotation[gc_index] & seurat@meta.data$tissue == "Salivary gland"] <- "SG GC"

  seurat@meta.data$grouped_annotation[seurat@meta.data$annotation %in% annotation[nbc_index] & seurat@meta.data$tissue == "Salivary gland"] <- "SG NBC"


  #Blood MBC, PC, GC
  seurat@meta.data$grouped_annotation[seurat@meta.data$annotation %in%  annotation[mbc_index] & seurat@meta.data$tissue == "Blood"] <- "PBMC MBC"

  seurat@meta.data$grouped_annotation[seurat@meta.data$annotation %in% annotation[plasma_index] & seurat@meta.data$tissue == "Blood"] <- "PBMC Plasma"

  seurat@meta.data$grouped_annotation[seurat@meta.data$annotation %in% annotation[gc_index] & seurat@meta.data$tissue == "Blood"] <- "PBMC GC"

  seurat@meta.data$grouped_annotation[seurat@meta.data$annotation %in% annotation[nbc_index] & seurat@meta.data$tissue == "Blood"] <- "PBMC NBC"

  #Other, non-B subsets
  seurat@meta.data$grouped_annotation[is.na(seurat@meta.data$grouped_annotation)] <- "Other"

  return(seurat)
}

#' A function to make volcano plots based on multiple marker dfs
#'
#' This function produces a patchwork of volcano plots based on a list of markers found by running the Seurat FindMarkers() function multiple times.
#' For example, this function can be used to plot multiple volcano plots side-by-side as a result of comparing two groups of cells within a cluster, for all clusters.
#' The thresholds for the adjusted p value is set to 0.05, the threshold for the avg_log2FC is set to 0.25.
#' @param marker_list A named list of the data.frame outputs of the Seurat FindMarkers(), FindAllMarkers() or FindConservedMarkers() function.
#' @param display_n An integer stating the amount of (significant) gene names that should be displayed on each side of the volcano plot.
#' @keywords seurat, volcano, marker, ggplot
#' @returns A list of ggplot objects with the volcano plots, as well as a figure displaying them side-by-side.
#' @export
#' @examples markers_1 <- FindMarkers(ident.1, ident.2)
#' markers_2 <- FindMarkers(ident.3, ident.4)
#' marker_list <- list(markers_1, markers_2)
#' plot <- volcano(marker_list, display_n = 10)
#'
volcano <- function(marker_list, display_n){
  stopifnot(require(cowplot))
  stopifnot(require(ggplot2))

  plot_list <- list()

  for (i in seq_along(marker_list)){

    #Label assignment
    marker_list[[i]]$sig <- "no"
    marker_list[[i]]$sig[marker_list[[i]]$p_val_adj <= 0.05] <- "yes"

    marker_list[[i]]$label <- NA
    marker_list[[i]] <- marker_list[[i]][order(marker_list[[i]]$avg_log2FC),]
    top_genes <- rbind(head(marker_list[[i]], display_n), tail(marker_list[[i]], display_n))
    top_genes <- top_genes[top_genes$sig == "yes",]

    matched <- which(rownames(marker_list[[i]]) %in% rownames(top_genes))
    marker_list[[i]]$label[matched] <- rownames(marker_list[[i]])[matched]

    marker_list[[i]] <- marker_list[[i]][order(marker_list[[i]]$p_val_adj),]
    top_genes <- head(marker_list[[i]], display_n)
    top_genes <- top_genes[top_genes$sig == "yes",]

    matched <- which(rownames(marker_list[[i]]) %in% rownames(top_genes))
    marker_list[[i]]$label[matched] <- rownames(marker_list[[i]])[matched]

    marker_list[[i]]$change <- "none"
    marker_list[[i]]$change[marker_list[[i]]$sig == "yes" & marker_list[[i]]$avg_log2FC > 0.25] <- "Up"
    marker_list[[i]]$change[marker_list[[i]]$sig == "yes" & marker_list[[i]]$avg_log2FC < -0.25] <- "Down"

    #Plotting
    ggplotcolors <- c("blue", "red", "gray")
    names(ggplotcolors) <- c("Down", "Up", "none")

    marker_list[[i]]$pct.1 <- as.numeric(marker_list[[i]]$pct.1)
    marker_list[[i]]$pct.2 <- as.numeric(marker_list[[i]]$pct.2)

    marker_list[[i]]$pct.diff <- abs(marker_list[[i]]$pct.1 - marker_list[[i]]$pct.2)

    #Generate a volcano plot indicating significant genes
    volcano <- ggplot(marker_list[[i]], aes(x=avg_log2FC, y=-log10(p_val_adj), col=change, label=label, size = pct.diff)) + geom_point() + ggtitle(names(marker_list)[[i]]) + geom_vline(xintercept=c(-0.25, 0.25), col="black") +
      geom_hline(yintercept=-log10(0.05), col="black") + scale_colour_manual(values = ggplotcolors) + geom_text_repel(box.padding=0.9)

    plot_list[[length(plot_list)+1]] <- volcano
  }
  suppressWarnings(cowplot::plot_grid(plotlist = plot_list))
  return(plot_list)
}

#' Make a UMAP of a specific clone
#'
#' This function produces a UMAP on which a specific clonotype is highlighted.
#' @param seurat A seurat object which has already been UMAP reduced and clustered, and has a metadata column indicating the clonotype of each cell
#' @param clone_number A character (vector) indicating the clonotype that needs to be highlighted on the UMAP
#' @keywords seurat, dimplot, umap, clonotype
#' @returns A UMAP with the indicated clonotypes highlighted.
#' @export
#' @examples plot <- plot_clone(seurat, c("1", "2"))
#'
plot_clone <- function(seurat, clone_number){
  stopifnot(require(Seurat))
  clone_number <- as.character(clone_number)

  if(!("clonotype_new" %in% colnames(seurat@meta.data))){
    if (grepl("clonotype", colnames(seurat@meta.data))){
      seurat@meta.data$clonotype_new <- seurat@meta.data[,which(grepl("clonotype", colnames(seurat@meta.data)))[[1]]]
    } else{
      stop("No clonotyping column was found.")
    }
  }

  if(!(clone_number %in% seurat@meta.data$clonotype_new)){
    stop("That clonotype number was not found in the Seurat object.")
  }

  seurat@meta.data$clone <- "NA"
  seurat@meta.data$clone[seurat@meta.data$clonotype_new %in% clone_number] <- clone_number
  clone_plot <- DimPlot(seurat, reduction = 'umap', label = F, repel = TRUE, label.size = 2.5, pt.size = 1,  group.by = "clone",
                        order = c(as.character(clone_number), "NA"), cols = c("grey", "red"))
  clone_plot
  return(clone_plot)
}

#' Filter out cells with abberant IGH chains
#'
#' This function filters out cells with two IGH chains or no IGH chain from an AIRR-formatted data frame.
#' @param data An AIRR-formatted data frame containing light chain and heavy chain entries for each cell, describing the genes they use, their CDR3, etc.
#' @keywords filtering, vdj, immcantation
#' @returns The same dataframe, but now with the mentioned cells filtered out. It also produces output on how many cells are removed.
#' @export
#' @examples data <- filtering_AIRR(data)
#'
filtering_AIRR <- function(data){
  data <- data[data$productive == TRUE,]

  # remove cells with multiple heavy chain
  multi_heavy <- table(filter(data, locus == "IGH")$cell_id)
  multi_heavy_cells <- names(multi_heavy)[multi_heavy > 1]

  print(paste0(length(multi_heavy_cells), " cells with multiple heavy chains have been removed"))

  data <- filter(data, !cell_id %in% multi_heavy_cells)

  # split cells by heavy and light chains
  heavy_cells <- filter(data, locus == "IGH")$cell_id
  light_cells <- filter(data, locus == "IGK" | locus == "IGL")$cell_id
  no_heavy_cells <- light_cells[which(!light_cells %in% heavy_cells)]

  print(paste0(length(no_heavy_cells), " cells with no heavy chains have been removed"))


  data <- filter(data, !cell_id %in% no_heavy_cells)
  return(data)
}

#' Extract all barcode labels for each tree in a list
#'
#' This function makes a list of all cell barcodes that are present in a phylogenetic tree.
#' @param db A phylogenetic tree in format "phylo" that contain trees of multiple cellular barcodes.
#' @keywords phylogenetic tree, barcodes
#' @returns A list, named with the same numbers as the tree numbers, of cell barcode vectors.
#' @export
#' @examples db <- readIghyml(".....igphyml-pass.tab", format = "phylo")
#' db_list <- extract_tree_sequences(db)
#'
extract_tree_sequences <- function(db){
  stopifnot(require(ape))

  output_list <- list()
  for (i in seq_along(db[["trees"]])){
    sequences <- db[["trees"]][[i]]$tip.label

    sequences <- sequences[!grepl("GERM", sequences)] #remove germline sequences

    output_list[[length(output_list)+1]] <- sequences
    names(output_list)[[i]] <- names(db[["trees"]])[[i]]
  }
  return(output_list)
}

#' Find autoreactive cells within a list of phylogenetic trees
#'
#' This function identifies all trees containing autoreactive cells, provided that a data.frame of autoreactive cell barcodes is given.
#' It also annotates every tree based on the metadata of the corresponding Seurat object.
#' @param db_seq A list of cellular barcodes in every tree in a list, as can be generated by the extract_tree_sequences() function.
#' @param donor A character vector indicating the donor of which the autoreactive cell barcodes should be used. Options are D001, D002, D003, D004LN, D004PG, D007 and AIM.
#' @param data_merged A data.frame containing all cell barcodes (of every donor) and their corresponding contig data.
#' @param autoreactive_cells A data.frame stating the barcode of each autoreactive cell and their donor of origin
#' @param seurat A seurat object containing the cluster annotation of each cell in the metadata column "annotation".
#' @keywords phylogenetic tree, autoreactivity
#' @returns A list of autoreactive cell barcodes, their cluster annotation, and the tree number they were found in.
#' @export
#' @examples db_seq <- extract_tree_sequences(db)
#' clones_1 <- find_autoreactive_clones(db_seq, donor = "D001", data_merged, autoreactive_cells, seurat)
#'
find_autoreactive_clones <- function(db_seq, donor, data_merged, autoreactive_cells, seurat){
  results <- list()
  data_merged$barcode <- paste0(data_merged$ident, "_", gsub("-1", "", data_merged$cell_id)) #make a barcode for each cell from the donor and sequence id

  if (!("annotation" %in% colnames(seurat@meta.data))){
    stop("Annotation column is missing.")
  }

  cluster_ident <- data.frame(barcode = rownames(seurat@meta.data), cluster = seurat@meta.data$annotation)

  #subset the AIRR object based on the donor
  if (donor == "D001"){
    data <- data_merged[data_merged$ident == "D001",]
  } else if (donor == "D002"){
    data <- data_merged[data_merged$ident == "D002",]
  } else if (donor == "D003"){
    data <- data_merged[data_merged$ident == "D003",]
  } else if (donor == "D004LN"){
    data <- data_merged[data_merged$ident == "D004LN",]
  } else if (donor == "D004PG"){
    data <- data_merged[data_merged$ident == "D004PG",]
  } else if (donor == "D007") {
    data <- data_merged[data_merged$ident == "D007",]
  } else if (donor == "B005-1") {
    data <- data_merged[data_merged$ident == "B005-1",]
  } else if (donor == "B005-2") {
    data <- data_merged[data_merged$ident == "B005-2",]
  } else if (donor == "B007-1") {
    data <- data_merged[data_merged$ident == "B007-1",]
  } else if (donor == "B007-2") {
    data <- data_merged[data_merged$ident == "B007-2",]
  } else {
    print("No entries found for that donor.")
  }


  for (i in seq_along(db_seq)){ #for each tree
    check_list <- list()
    sequences <- db_seq[[i]] #get the sequences out
    indexes <- which(data$sequence_id %in% sequences) #get the index numbers out of the sequences in the AIRR object

    barcodes <- data$barcode[indexes] #grab the corresponding barcodes

    result_pertree <- data.table

    for (j in barcodes){ #for each barcode
      index_bar <- which(barcodes %in% j)
      seq_index <- indexes[[index_bar]]

      check <- autoreactive_cells[autoreactive_cells$barcode %in% j,] #check if it occurs in the list of autoreactive cells

      if (nrow(check) > 0){
        check$sequence_id <- data$sequence_id[seq_index]
        check_list[[length(check_list)+1]] <- check #if the barcode does appear in the list of autoreactive cells, add it to the output

      } else{
        next #if the barcode does not appear, go to the next sequence
      }

    }
    output <- rbindlist(check_list)

    if (nrow(output) == 0) { #if none of the sequences are autoreactive, skip the tree
      next
    }

    #take cluster identities from seurat object
    output <- left_join(output, cluster_ident, by = "barcode")

    results[[length(results)+1]] <- output
    names(results)[[length(results)]] <- names(db_seq)[[i]]
  }



  return(results)
}

#' Generate a phylogenetic tree plot
#'
#' This function generates a plot of a phylogenetic tree, indicating any possible autoreactive sequences in the tree, and also adding the cluster annotation from Seurat.
#' For each donor, a corresponding db object and autoreactive (auto) data frame should exist. For example, for D001, a phylo tree named db1 should exist, as well as auto1, which is the result of running the find_autoreactive_clones() function for D001.
#' @param tree_number A character of the tree number identifier (as given in the names of the tree list) that is to be plotted.
#' @param seurat A seurat object with the cluster annotation of each cell in a metadata column called "annotation".
#' @param ident The donor ident to be used when determining if the sequences are autoreactive or not.
#' @keywords phylogenetic tree, autoreactivity, plot
#' @returns A phylogenetic tree plot of the indicated tree. The cluster identities for each cell are plotted as a bar to the side. A legend for this cluster annotation is added in the upper left corner.
#' A cluster identity of NA means that the cell was not found back in the suerat object.
#' The cells are coloured based on their autoreactivity: black if non-autoreactive, green in RF, pink if Ro60, gold if Ro52, purple if both Ro52 and Ro60, red if La, and cyan if the cell is autoreactive but does not match any of these combinations.
#' @export
#' @examples plot <- phylo_tree("123_1", ident = "D003", seurat)
#'
phylo_tree <- function(tree_number, seurat = NULL, ident, annotation = T){

  if (annotation == T & is.null(seurat)){
    stop("A Seurat object needs to be defined for tree annotation. If you do not need cluster annotations in your tree, then set annotation to FALSE.")
  }

  grid.newpage()


  if (ident == "D001"){
    tree_object <- db1
    auto <- auto_1
  } else if (ident == "D002"){
    tree_object <- db2
    auto <- auto_2
  } else if (ident == "D003"){
    tree_object <- db3
    auto <- auto_3
  } else if (ident == "D004LN"){
    tree_object <- db4ln
    auto <- auto_4ln
  } else if (ident == "D004PG"){
    tree_object <- db4pg
    auto <- auto_4pg
  } else if (ident == "D007"){
    tree_object <- db7
    auto <- auto_7
  } else if (ident == "B005-1"){
    tree_object <- db0051
    auto <- auto_b0051
  } else if (ident == "B005-2"){
    tree_object <- db0052
    auto <- auto_b0052
  } else if (ident == "B007-1"){
    tree_object <- db0071
    auto <- auto_b0071
  } else if (ident == "B007-2"){
    tree_object <- db0072
    auto <- auto_b0072
  } else{
    stop("Donor not found.")
  }

  if (!(tree_number %in% names(tree_object$trees))){
    stop("Tree number not found for this donor.")
  }

  tree <- tree_object$trees[[as.character(tree_number)]]
  tree <- as(tree, "phylo4")

  data <- as.data.frame(tree@label)
  colnames(data) <- c("sequence_id")

  if (tree_number %in% auto$tree){
    auto <- auto[order(match(auto$sequence_id, data$sequence_id)),]


    data <- left_join(data, auto, by = "sequence_id")
    data <- data[order(match(data$sequence_id, tree@label)),]

    col <- seq(1:nrow(data))

    for (i in 1:nrow(data)){
      if (!is.na(data[i,c("antigen")])){
        if (grepl("RF", data[i,c("antigen")])){
          col[i] <- "darkgreen"
        } else if (data[i,c("antigen")] == "Ro60"){
          col[i] <- "hotpink"
        } else if (data[i,c("antigen")] == "Ro52"){
          col[i] <- "gold"
        } else if (data[i,c("antigen")] == "La"){
          col[i] <- "red"
        } else if (data[i,c("antigen")] == "Ro60,Ro52"){
          col[i] <- "purple"
        } else {
          col[i] <- "cyan"
        }

      } else{
        col[i] <- "black"

      }
    }
  } else {
    col <- rep("black", times = nrow(data))
    seurat_clusters <- data.frame(cell_barcode = rownames(seurat@meta.data), cluster = seurat@meta.data$annotation)
    seurat_clusters <- seurat_clusters[grepl(ident, seurat_clusters$cell_barcode),]
    seurat_clusters$barcode <- gsub(paste0(ident, "_"), "", seurat_clusters$cell_barcode)

    data$barcode <- gsub("-1_contig_2|-1_contig_1", "", data$sequence_id)


    # Sort seurat_clusters based on the order
    seurat_clusters <- seurat_clusters[seurat_clusters$barcode %in% data$barcode,]
    data <- left_join(data, seurat_clusters, by = "barcode")
    data[,c("cell_barcode", "barcode")] <- NULL
  }

  for (i in 1:nrow(data)){
    for (j in 1:nrow(data)){
      barcode_1 <-  gsub("-1_contig_2|-1_contig_1|-1_contig_1_1|-1_contig_2_1", "", data$sequence_id[[i]])
      barcode_2 <-  gsub("-1_contig_2|-1_contig_1|-1_contig_1_1|-1_contig_2_1", "", data$sequence_id[[j]])

      if (barcode_1 == barcode_2){
        cols <- col[c(i,j)]
        cols <- cols[cols != "black"]

        if(length(unique(cols)) == 1){
          col[c(i,j)] <- cols
        }
      }
    }
  }


  #generate a coloured sidebar for the cluster
  if(annotation == T){
    colors <- distinctColorPalette(k = length(unique(seurat$annotation)))
    unique_clusters <- unique(seurat$annotation)
    cluster_col <- data.frame(cluster = unique_clusters, cluster_color = colors)
    cluster_col <- rbind(cluster_col, data.frame(cluster = NA, cluster_color = "#000000"))

    cluster_col <- cluster_col[cluster_col$cluster %in% data$cluster,]

    data <- inner_join(cluster_col, data, by = "cluster")
    data$cluster_color[grepl("GERM", data$sequence_id)] <- NA
    data <- data[order(match(data$sequence_id, tree@label)),]

    colors <- data$cluster_color
    colors <- colors[!is.na(colors)]
  }


  #Make the tree plot
  phylo <- phyloXXYY(tree)

  annot_bar <- grid::viewport(width = 0.8, height = 0.8)
  grid::pushViewport(annot_bar)

  # Add the colored rectangles
  if (annotation == T){
    sequence_y <- seq(-0.1, 1.02, by = 1/length(colors))
    sequence_y <- rev(sequence_y)



    for (i in 1:length(colors)) {
      grid.rect(x = 1.1, y = sequence_y[[i]], width = 0.05, height = 1/length(colors),
                default.units = "npc", just = "center", gp = gpar(fill = colors[i]))
    }


    #Add in the makeshift legend


    for (i in 1:nrow(cluster_col)){
      grid.rect(x = 0 , y = 0.9 + (i/30), width = 0.05, height = 0.04,
                default.units = "npc", just = "center", gp = gpar(fill = cluster_col$cluster_color[i]))

      grid.text(cluster_col$cluster[[i]], x = 0.1, y = 0.9 + (i/30),
                just = "centre", hjust = NULL,  default.units = "npc",
                gp = gpar())

    }
  }

  #Plot the tree
  plot <- plotOneTree(phylo, type = "phylogram", show.tip.label = TRUE,
                      show.node.label = FALSE, tip.color = col, edge.color = 'black',
                      edge.width = 1, node.color = "black", rot = 0)


  # Exit the viewport
  grid::popViewport()

  return(plot)

}



#Bugfixed some scRepertoire code
combineExpression.new <- function(df, sc, cloneCall="gene",
                                  chain = "both", group.by="none",
                                  proportion = TRUE, filterNA = FALSE,
                                  cloneTypes=c(Rare = 1e-4, Small = 0.001,
                                               Medium = 0.01, Large = 0.1, Hyperexpanded = 1),
                                  addLabel = FALSE) {
  options( dplyr.summarise.inform = FALSE )
  stopifnot(require(scRepertoire))
  cloneTypes <- c(None = 0, cloneTypes)
  df <- checkList(df)
  cloneCall <- theCall(cloneCall)
  Con.df <- NULL
  meta <- grabMeta(sc)
  cell.names <- rownames(meta)
  if (group.by == "none") {
    for (i in seq_along(df)) {
      if (chain != "both") {
        df[[i]] <- off.the.chain(df[[i]], chain, cloneCall)
      }
      data <- data.frame(df[[i]], stringsAsFactors = FALSE)
      data2 <- unique(data[,c("barcode", cloneCall)])
      data2 <- na.omit(data2[data2[,"barcode"] %in% cell.names,])
      if (proportion == TRUE) {
        data2 <- data2 %>% group_by(data2[,cloneCall]) %>%
          summarise(Frequency = n()/nrow(data2))
      } else {
        data2 <- data2 %>% group_by(data2[,cloneCall]) %>%
          summarise(Frequency = n())
      }
      colnames(data2)[1] <- cloneCall
      data <- merge(data, data2, by = cloneCall, all = TRUE)
      data <- data[,c("barcode", "CTgene", "CTnt",
                      "CTaa", "CTstrict", "Frequency")]
      Con.df <- rbind.data.frame(Con.df, data)
    }
  } else if (group.by != "none") {
    data <- data.frame(bind_rows(df), stringsAsFactors = FALSE)
    data2 <- na.omit(unique(data[,c("barcode", cloneCall, group.by)]))
    data2 <- data2[data2[,"barcode"] %in% cell.names, ]
    data2 <- as.data.frame(data2 %>% group_by(data2[,cloneCall],
                                              data2[,group.by]) %>% summarise(Frequency = n()))
    colnames(data2)[c(1,2)] <- c(cloneCall, group.by)
    x <- unique(data[,group.by])
    for (i in seq_along(x)) {
      sub1 <- subset(data, data[,group.by] == x[i])
      sub2 <- subset(data2, data2[,group.by] == x[i])
      merge <- merge(sub1, sub2, by=cloneCall)
      if (proportion == TRUE) {
        merge$Frequency <- merge$Frequency/length(merge$Frequency)
      }
      Con.df <- rbind.data.frame(Con.df, merge)
    }
    nsize <- Con.df %>% group_by(Con.df[,paste0(group.by, ".x")])  %>% summarise(n = n())
  }

  Con.df$cloneType <- NA
  for (x in seq_along(cloneTypes)) { names(cloneTypes)[x] <-
    paste0(names(cloneTypes[x]), ' (', cloneTypes[x-1],
           ' < X <= ', cloneTypes[x], ')') }
  for (i in 2:length(cloneTypes)) { Con.df$cloneType <-
    ifelse(Con.df$Frequency > cloneTypes[i-1] & Con.df$Frequency
           <= cloneTypes[i], names(cloneTypes[i]), Con.df$cloneType) }
  PreMeta <- unique(Con.df[,c("barcode", "CTgene", "CTnt",
                              "CTaa", "CTstrict", "Frequency", "cloneType")])
  dup <- PreMeta$barcode[which(duplicated(PreMeta$barcode))]
  PreMeta <- PreMeta[PreMeta$barcode %!in% dup,]
  rownames(PreMeta) <- PreMeta$barcode
  if (group.by != "none" && addLabel) {
    location <- which(colnames(PreMeta) == "Frequency")
    colnames(PreMeta)[location] <- paste0("Frequency.", group.by)
  }
  if (inherits(x=sc, what ="Seurat")) {
    if (length(which(rownames(PreMeta) %in%
                     rownames(sc[[]])))/length(rownames(sc[[]])) < 0.01) {
      warning("< 1% of barcodes match: Ensure the barcodes in
            the Seurat object match the
            barcodes in the combined immune receptor list from
            scRepertoire - most common issue is the addition of the
            prefixes corresponding to `samples` and 'ID' in the combineTCR/BCR()
            functions")
    }
    col.name <- names(PreMeta) %||% colnames(PreMeta)
    sc[[col.name]] <- PreMeta
  } else {
    rownames <- rownames(colData(sc))
    if (length(which(rownames(PreMeta) %in%
                     rownames))/length(rownames) < 0.01) {
      warning("< 1% of barcodes match: Ensure the barcodes
          in the SingleCellExperiment object match the
          barcodes in the combined immune receptor list from
          scRepertoire - most common issue is the addition of the
          prefixes corresponding to `samples` and 'ID' in the combineTCR/BCR()
          functions") }
    colData(sc) <- cbind(colData(sc), PreMeta[rownames,])[, union(colnames(colData(sc)),  colnames(PreMeta))]
    rownames(colData(sc)) <- rownames
  }
  if (filterNA == TRUE) { sc <- filteringNA(sc) }
  return(sc)
}


"%!in%" <- Negate("%in%")

#Use to shuffle between chains
off.the.chain <- function(dat, chain, cloneCall) {
  chain1 <- toupper(chain) #to just make it easier
  if (chain1 %in% c("TRA", "TRG", "IGH")) {
    x <- 1
  } else if (chain1 %in% c("TRB", "TRD", "IGL")) {
    x <- 2
  } else {
    warning("It looks like ", chain, " does not match the available options for `chain = `")
  }
  dat[,cloneCall] <- str_split(dat[,cloneCall], "_", simplify = TRUE)[,x]
  return(dat)
}


#Ensure df is in list format
checkList <- function(df) {
  df <- if(is(df)[1] != "list") list(df) else df
  return(df)
}



#This is to grab the meta data from a seurat or SCE object
#' @importFrom SingleCellExperiment colData
#' @importFrom SeuratObject Idents
grabMeta <- function(sc) {
  if (inherits(x=sc, what ="Seurat")) {
    meta <- data.frame(sc[[]], slot(sc, "active.ident"))
    colnames(meta)[length(meta)] <- "ident"
  } else if (inherits(x=sc, what ="SummarizedExperiment")){
    meta <- data.frame(colData(sc))
    rownames(meta) <- sc@colData@rownames
    clu <- which(colnames(meta) == "ident")
    colnames(meta)[clu] <- "ident"
  }
  return(meta)
}

#Filtering NA contigs out of single-cell expression object
#' @import dplyr
#' @importFrom SingleCellExperiment colData
filteringNA <- function(sc) {
  meta <- grabMeta(sc)
  evalNA <- data.frame(meta[,"cloneType"])
  colnames(evalNA) <- "indicator"
  evalNA <- evalNA %>%
    transmute(indicator = ifelse(is.na(indicator), 0, 1))
  rownames(evalNA) <- rownames(meta)
  if (inherits(x=sc, what ="cell_data_set")){
    colData(sc)[["evalNA"]]<-evalNA
    return(sc[, !is.na(sc$cloneType)])
  }else{
    col.name <- names(evalNA) %||% colnames(evalNA)
    sc[[col.name]] <- evalNA
    sc <- subset(sc, cloneType != 0)
    return(sc)
  }
}


# This is to help sort the type of clonotype data to use
theCall <- function(x) {
  if (x %in% c("CTnt", "CTgene", "CTaa", "CTstrict")) {
    x <- x
  }else if (x == "gene" | x == "genes") {
    x <- "CTgene"
  } else if(x == "nt" | x == "nucleotide") {
    x <- "CTnt"
  } else if (x == "aa" | x == "amino") {
    x <- "CTaa"
  } else if (x == "gene+nt" | x == "strict") {
    x <- "CTstrict"
  }
  return(x)
}



