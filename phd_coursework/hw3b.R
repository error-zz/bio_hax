
# problem 3.3.b

f_b_string <- function (b_string) {
  returner <- 0
  b_string <- trimws(b_string)
  b_string_chars <- strsplit(b_string, " ")[[1]]
  counter <- 1
  for (char in b_string_chars) {
    char <- as.numeric(char)
    returner <- returner + 2^(counter*char)
    string <- paste("C ",char," : ",returner)
    # print(string)
    counter <- counter + 1
  }
  string <- paste("S ",b_string," : ",returner)
  # print(string)
  return(returner)
}

new_f_z <- function (b_string){
  string <- paste("Evaluating ",b_string)
  # print(string)
  f_of_b <- f_b_string(b_string)
  returner <- (1-a)/(1+a)*(a^abs(z-f_of_b))
  string <- paste("R ",returner)
  # print(string)
  return (returner)
} 

# hard-coded (problem formulation)
regression_count <- 50 # modify for each regression count test
z <- 128
a <- 0.25
b_max <- 10

# main
num_sum <- 0
den_sum <- 0
b_sum <- 0
all_sum <- 0
den_ct <- 0

# for each regression
test <- array(0,dim=c(1,regression_count))
for (r in 1:regression_count){
  #populate simulated bits
  b_string <- ""
  den_sum <- 0
  num_sum <- 0
  compare_char <- ""
  for (b in 1:b_max){
    new_char <- sample(0:1, 1)
    b_string <- paste (b_string, new_char)
    if (b == 7){ compare_char <- new_char }
  }
  calc <- new_f_z(b_string)
  string <- paste("CC: ", compare_char)
  compare_char <- as.numeric(compare_char)
  # print(string)
  if (compare_char == TRUE){ # indicator variable
      # print("!")
      num_sum <- num_sum + (calc*0.5)
  }
  den_sum <- den_sum + calc
  tmp <- num_sum/den_sum
  if (is.nan(tmp)){ test[r] <- 0.5 } # clean up infinity cases
  else{ test[r] <- tmp }
  all_sum <- all_sum + test[r]
}
# string <- paste("All test cases: ", test)
# print(string)
avg <- all_sum/regression_count
string <- paste("Average of all tests: ", avg)
print(string)
# barplot(test)