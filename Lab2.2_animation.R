#*******************#
# Lab2: Annimations
#*******************#

#===============================#
# define a plot function of time
#===============================#

# the funtion that plots the spatial map as a function of time

Tmax_t <- function(tau) {
  Tmax_sub <- filter(Tmax, t == tau)
  
  ggplot(Tmax_sub) +
    geom_point(aes(x = lon, y = lat, colour = z),
               size = 4) +
    col_scale(name = "z", limits = c(40, 110)) + 
    theme_bw()
}


#==========================#
# plot the seq of animation
#==========================#

# Next, construct a function that plots
# the data for every day in the data sets

# the function generates the animation within the HTML webpage

saveHTML(expr = # contains the all the images to plot)


gen_anim <- function() {
  for(t in lim_t[1]:lim_t[2]) {
    plot(Tmax_t(t))
    
  }
}           

ani.options(interval = 0.2) # 0.2s between frames
saveHTML(gen_anim(),
         autoplay = FALSE,
         loop = FALSE,
         verbose = FALSE,
         outdir = ".",  # save to current dir
         single.opts = " 'controls': ['first', 'previous', 
         'play', 'next', 'last', 'loop',
         'speed'], 'delayMin' : 0", 
         htmlfile = "NOAA_anim.html")

# teh max T clearly drifts from west to east 
# at several points during the animation

# suggest a dynamic spatio-temp model could
# capture this drift could provide a good fit
# to these data
























