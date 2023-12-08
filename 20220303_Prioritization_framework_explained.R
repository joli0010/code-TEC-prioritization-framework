#Load library:
library(tidyverse)

##### This is the code for Scenario 2a: Combined-approach

#setting up the data:
setwd()

Cost_TAS <- read.csv("20230418_FINAL_COSTS_low_best_upper.csv") #dataframe with costs per threat abatement strategy
areas <- read.csv("20230418_areas1476_final.csv") #area of the distribution of each species and TECs (to manage), in km2
sptec_int <- read.csv("20230418_sp_to_TECspisect_1476_final.csv") #overlap between each species and TEC. Each row represents an overlap, with the area of overlap and the proportion of the overlap relative to the total distribution of the species and the TEC

all_threats_matrix_or <- read.csv("20230418_MATRIX_2b_1476_final.csv") #TEC/sp - threat abatement strategy matrix
all_threats_matrix_or$id_sp <- as.character(all_threats_matrix_or$id_sp)
all_threats_matrix_or <- all_threats_matrix_or %>% mutate(weight = 1, .before = 4)


# Calculating all species benefits, considered individually
U <- matrix(1,nrow=nrow(all_threats_matrix),ncol=nrow(Cost_TAS)) 

benefits_scenario2a <- function(){
  rs <- data.frame(species_id = all_threats_matrix$id_sp, B = NA, C = NA) #RESULTS/OUTPUT; B= BENEFIT; C= COST
  for (i in 1:nrow(rs)){
    overlaps <- subset(sptec_int, sp_a == all_threats_matrix$id_sp[i]) #subsets the species I'm interested in
    Tik = as.numeric(all_threats_matrix[i,5:18]) #it tells Tik where the threats are positioned
    Uik = U[i,] #creates a row with x columns and only 1 values
    if (nrow(overlaps)==0){ #only dealing with species that have NO overlaps 
      cat("No overlaps found for species ", i, ", ",all_threats_matrix$id_sp[i], "\n") 
      rs[i,"B"] =  all_threats_matrix$weight[i] * sum(Tik * Uik) / sum(Tik) #benefit f(x) without overlap -benefit 1 just to self
      rs[i,"C"] = sum(Tik * Cost_TAS$Best_cost_scenario) * subset(areas, id_sp == all_threats_matrix$id_sp[i])$KM2[1] #cost function (same as the one below)
      next
    }
    
    mtch <- match(overlaps$sp_b, all_threats_matrix$id_sp) #shows the order in the df of the species that overlapped
    overlaps$W <- all_threats_matrix$weight[mtch] 
    Tjk = as.matrix(all_threats_matrix[mtch,5:18])  #creates a matrix for the binary variable with the species that overlapped and threat positions (columns represent threats)
    overlaps$total_threats <- apply(Tjk,1,sum) #total number of threats, by summing the threats in Tjk
    Ujk = U[mtch,] # Ujk is a matrix with threats in columns and species that overlapped in rows - only 1s inside
    overlaps$common_threats <- as.numeric((Ujk*Tjk) %*% Tik) #this multiplier obtains whether threats are in common or not
    overlaps$benefit_to_j <-  with(overlaps, prop_overlap * W * common_threats / total_threats) #this is essentially the benefit function
    overlaps <- subset(overlaps,!is.na(W) & total_threats>0) ##some species have total_threats=0 given we dont have costs for all threats, so probably some species have other threats that are not in the costs list
    
    rs[i,"B"] = sum(overlaps$benefit_to_j) ##total benefit
    rs[i,"C"] = sum(Tik * Uik * Cost_TAS$Best_cost_scenario) * subset(areas, id_sp == all_threats_matrix$id_sp[i])$KM2[1]
    rs[i,"B"] =  rs[i,"B"] + all_threats_matrix$weight[i] * sum(Tik * Uik) / sum(Tik) ## to stop double counting benefit of self
  }     ## automatic benefit to initial sp, same as the benefit inside the if statement, only this is for species that do overlap with others
  rs$E = rs$B/rs$C
  return(rs)
}


# Run code
combined <- benefits_scenario2a()
write_csv(combined, "combined.csv")


## To remove top species/TECs for reaching budget:
adjust_U <- function(U,species){ 
  i = which(all_threats_matrix$id_sp == species) # for species located in the threat matrix which the same id-name, call it i
  U[i,] <- 0  #then for i (determined above), transform to 0 the values in U matrix
  overlaps <- subset(sptec_int, sp_a == all_threats_matrix$id_sp[i]) #sub-setting the overlap table to only the species i
  mtch <- match(overlaps$sp_b, all_threats_matrix$id_sp) #shows vector position that matched that condition
  overlaps <- overlaps[!is.na(mtch),] #deleting NA values
  mtch <- mtch[!is.na(mtch)]
  Tik = as.numeric(all_threats_matrix[i,5:18]) #threats applicable to sp i
  U[mtch,] <- U[mtch,] * (1-outer(overlaps$prop_overlap,Tik,"*")) #overlaps$prop_overlap,Tik,"*" ==> for the proportion of overlaps,
                                  #multiply each value by Tik, which is the variable that tells you whether a column is a threat for species i or not
                                  #it creates a new matrix, only for the species it matched (overlapped) with the values of the proportion of overlaps under the columns of the threats for species i,
                                  #so it does 1 minus the proportion of overlap of these species (IT IS NOT THE PROPORTION OF THREATS), and U[mtch,]<- what it does is to add the 1-outer... to the whole U matrix, to then run the code again
  return(U)
}

#### When I do the next round of removing the second top species, now the values for this second species get added onto the U matrix
#(values=proportion of overlap of the overlapping species in the columns of the threats for species i)
#if there is an overlapping species (species j) that overlaps with BOTH top species (species i), then the proportion of overlap of the second top species is reduced from the value already there. 
#So e.g., if a koala and bilby overlapped by 60% of the bilby distribution, then 0.4 (1-0.6=0.4) would be added in matrix U under the columns of the threats for koala;
#if then the second umbrella species i, a parrot, also overlapped with the bilby, let's say by 20%, then the U matrix would change this under the same columns as above and now it would look like 0.4-0.2 = 0.2 BECAUSE now that's the remaining proportion of overlap that is NOT being managed by the koala and the parrot (top species 1 and 2)
#PLUS-- in case the second top species (parrot) shared other threats with the bilby (not included by the koala), then these columns would have the 1-0.2 now, so 0.8 under the new threats that are considered by managing the parrot.
#So then every time we exclude a top species because we're already considering its management, the proportion of overlap and the threats that are being managed for the benefiting species are being accounted for, and will not be double counted.

#WHAT HAPPENS WITH THE COSTS:
#When the costs get recalculated after removing a top species, Uik has changed, as instead of 1s, it might have a proportion of overlap of what has not been managed yet (not only for the top species but also for all other species that are considered as species i in the iteration)
#This means that now the code is costing ONLY for the proportion that has NOT been managed yet And when the proportion of species j that has been managed by the top species (now removed), is 1 (thus, the TOTAL area of species j), that's when the costs for species j are now 0 and CE is Inf, because all threats have been managed for it's total distribution
    

#FUNCTION TO REMOVE THREATS FOR TOP SPECIES THAT HAVE ALREADY BEEN MANAGED:
adjust_matrix <- function(all_threats_matrix,species){ 
  i = which(all_threats_matrix$id_sp == species) # for species located in the threat matrix which the same id-name, call it i
  all_threats_matrix[i, c(4:18)] <- 0  #then for i (determined above), transform to 0 the threat values in threats_matrix
  return(all_threats_matrix)
}


#TO REMOVE TOP SPECIES WITH A WHILE LOOP:
priority_combined <- data.frame()
scenario2a_bestcost_E <- combined
budget <- 125660000
total_cost <- 0

while(total_cost <= budget){
  
  # 1) order species/TECs by cost-effectiveness - descending order
  scenario2a_bestcost_E <- scenario2a_bestcost_E %>% arrange(-E)
  
  # 2) Removing Inf values in E
  scenario2a_bestcost_E <- scenario2a_bestcost_E %>% filter(!E == Inf)
  
  # 3) Add top species to priority_species df
  priority_combined <- rbind(priority_combined, scenario2a_bestcost_E[1, 1:4])
  
  # 4) Add cost to total cost
  total_cost <- total_cost + scenario2a_bestcost_E[1, "C"]
  
  # 5) Remove top species
  U <- adjust_U(U, scenario2a_bestcost_E[1,1]) 
  
  # 6) Remove top species from the iteration - Tik should be 0, as all threats for top species have been managed - so changing matrix to 0
  all_threats_matrix <- adjust_matrix(all_threats_matrix, scenario2a_bestcost_E[1,1])
  
  # 7) Run function again, now without top species
  scenario2a_bestcost_E <- SPTECspecies_benefits()
}

write_csv(priority_combined, "priority_list_combined.csv") #this is the output file with the list of priority species and TECs. It includes surrogate and additional features.



