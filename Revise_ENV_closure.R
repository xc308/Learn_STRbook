#************#
# Advanced R
# chp 8 
#************#

search()# lists all parents of the global env. 
        # called search path as objs in these envs can be found at the top-level interactive workspace (globalenv())
        # contains one env for each attached package and objects attached
        # also contains special env called Atuoloads, used to save memory, load big datasets when needed. 
        # each new package is loaded with library(), it's inserted in betwn the global env and previous package at the top of search list

search()


#==================#
# Create a new env
#==================#

e <- new.env()
# defautl parent of new.env() is the env from which it's called, so global env here
parent.env(e)
# <environment: R_GlobalEnv>


#------------------------------#
# modify binding within this env
#------------------------------#

# treat it like a list
e$a <- 1
e$b <- 2

ls(e) # [1] "a" "b"

e$a # [1] 1

## ls() only list the names don't befin with .
  # need all.names = TRUE to show all names

e$.a <- 2
ls(e) # [1] "a" "b"
ls(e, all.names = T) # [1] ".a" "a"  "b" 


## Another way to view an env is ls.str(), 
  # better than ls()
  # as it shows each obj in the env

str(e) # <environment: 0x7f88e4492270> 
ls.str(e) 
# a :  num 1
# b :  num 2

ls.str(e, all.names = T)
# .a :  num 2
# a :  num 1
# b :  num 2


#----------------------------#
# Given a name, extract value
#----------------------------#

# $, [[ ]]: only look in one env
# get() : use scoping rules

e$c # NULL
e$c <- 3

e[["c"]] # [1] 3

get("c", envir = e) # [1] 3


#-----------------------#
# Deleting objs from env
#-----------------------#

# use rm() to remove the binding
rm("a", envir = e)

ls.str(e, all.names = T) 
# .a :  num 2
# b :  num 2
# c :  num 3


#---------------------------#
# if a binding exists in env
#---------------------------#

# default behavior: fowllow scoping rules to look in parent env
  # inherits = F, to stop this

x <- 10

exists("x", envir = e)  
# [1] TRUE
# x not in e, but in global, which is the parent of e
# so follow scoping rules


exists("x", envir = e, inherits = F)
# [1] FALSE


#-------------#
# compare env
#-------------#

# use indentical()
identical(globalenv(), environment()) # [1] TRUE



#========================#
# 8.2 Recursing over envs
#========================#

# Env forms a tree, so it's convenient to write
# recursive function

# understand pryr::where()

install.packages("pryr")
library("pryr")

# Given a name, where() finds the env where that name is defined

x <- 5
where("x") # <environment: R_GlobalEnv>

where("mean") # <environment: package:Matrix>

# has two args: 
  # name (str) to look for
  # env in which to start the search


where <- function(name, env = parent.env()) {
  if(identical(env, emptyenv())){
    # Base case, reached empty env
    stop("Can't find", name, call. = T)
  } else if (exists(name, envir = env, inherits = F)) {
    # success case
    env
  } else {
    # recursive case, try parent env
    where(name, parent.env())
  }
}


x <- 0 
f <- function(){
  x <<- 1
}

f()
x


new_counter <- function(){
  i <- 0
  function(){
    i <<- i + 1 
    i
  }
}

i <- 0
new_counter2 <- function() {
  i <<- i + 1
  i 
}
new_counter2()
# 1
# 2
# 3
#...

new_counter3 <- function() {
  i <- 0
  funtion() {
    i <- i + 1
    i
  }
}


i <- 0
new_counter4 <- function() {
  i <- i + 1
  i
}

new_counter4()









































































