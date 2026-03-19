# Calculates Jaccard similarity between two users (columns) in a binary matrix
#Explanation of Inputs:
  
# extract gene expression and change to binary; expression value greater than 0 = 1; 

#`M`: A binary matrix/data frame where rows represent items and columns represent users (or samples/features), with values of 0 or 1.
#`user1`: The column name (or index) for the first user/sample.
#`user2`: The column name (or index) for the second user/sample.



jaccard_2 <- function(M, user1, user2) {
  # Sum each row for the two specified users (columns)
  sums = rowSums(M[, c(user1, user2)])
  
  # Count the number of rows where both users have 1 (intersection)
  similarity = length(sums[sums == 2])
  # Count the number of rows where at least one user has 1 (union minus intersection), then add intersection
  total = length(sums[sums == 1]) + similarity
  
  # Return Jaccard index: intersection divided by union
  similarity/total
}

