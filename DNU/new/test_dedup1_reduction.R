index_tmp <- index_dedup1_general_use[[1]][[19]]

index_tmp2 <- index_tmp
for(i in 1:ncol(index_tmp)){
  it <- index_tmp[, i]
  for(j in 1:2){
    if(it[2 * j - 1] <= it[2 * j]){
      index_tmp2[2 * j - 1, i] <- it[2 * j - 1]
      index_tmp2[2 * j, i] <- it[2 * j]
    }else{
      index_tmp2[2 * j - 1, i] <- it[2 * j]
      index_tmp2[2 * j, i] <- it[2 * j - 1]
    }
  }
}
bool_duplicated_adopt <- !duplicated(t(index_tmp2))
index_dedup1_tmp <- index_tmp2[, bool_duplicated_adopt, drop = FALSE]

ncol(index_tmp)
ncol(index_dedup1_tmp)
