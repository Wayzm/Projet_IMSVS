# ----- IMSVS ----- #
library('plot.matrix')

## PARAMS

y_max = 50
x_max = 50
DThreshold = 50
IThreshold = 5
init_health = 20
init_food_value = 5
n_bacts = 100
r_bact = c(3,2,1) # ratio initial de bactéries B:C:D

lyso_time = 600 # lysosyme is released at steady-state
lyzo_dmg = 1 # zero pour desactiver

cost_of_living = 5 # health penalty per iteration
division_factor = 1 # health penalty on division (%percentage)

food_per_itt = 1 # number of food to add on each cell per iteration
food_delay = 1 # refill food delay

timestamp = 1000 # duration of simulation
plotfreq = 100
set.seed(4)

#autre param dans Comsumption 
# -> lactate production by B
# -> prevent or allow multiple feeding on same cell

steady_state_size = 50 # taille de l'intervale à moyenner pour obtenir Bs,Cs,Ds
color=c("#333333", "#BF75F9", "#35AEF7", "#CFCB3A")

##########
# FUNCTIONS

# Placing  bacteria
Init_bacteria<-function(n_bacts, r_bact)
{
  Health <<- matrix(data = 0, ncol =  x_max, nrow = y_max)
  Type <<- matrix(data = 0, ncol =  x_max, nrow = y_max)
  
  id = sample(1:(x_max*y_max), n_bacts, replace = FALSE)
  
  # setup nombre initial de bacterie de chaque type 
  r_bact = r_bact * n_bacts / sum(r_bact)
  r_bact[2] = r_bact[2] + r_bact[1]
  
  for (i in 1:n_bacts) {
    x = ((id[i]-1) %% y_max) + 1
    y = floor((id[i]-1) / y_max) + 1
    Health[x,y] <<- init_health
    
    if     (i <= r_bact[1]) Type[x,y] <<- 1 # 'B'
    else if(i <= r_bact[2]) Type[x,y] <<- 2 # 'C'
    else                    Type[x,y] <<- 3 # 'D'
  }
}

# Placing the food
Placing_food<-function()
{
  Food <<- array(rep(init_food_value, x_max * y_max * 6), dim=c(x_max, y_max, 6))
}

Add_food<-function()
{
  Food <<- Food + food_per_itt
}

# Cell Division
Cell_Division<-function(x, y)
{
  original_health = floor(Health[x, y] / 2 * division_factor)
  typeb<-Type[x, y]
  
  xs <- c(x - 1, x    , x + 1, x    , x - 1, x - 1, x + 1, x + 1)
  ys <- c(y    , y + 1, y    , y - 1, y - 1, y + 1, y + 1, y - 1)
  
  # Join edge of the grid
  for(i in 1:8){
    if(xs[i] < 1 || xs[i] > x_max)
      xs[i] = ((xs[i] - 1)%%x_max) + 1
    if(ys[i] < 1 || ys[i] > y_max)
      ys[i] = ((ys[i] - 1)%%y_max) + 1
  }
  
  # We get the number of empty cells and their coords
  empty_cells = 0 # Nombre de cellule vide 
  cells_coords <- matrix(data = 0, ncol =  8, nrow = 2)
  
  for(i in 1:8)
  {
    xx = xs[i]
    yy = ys[i]
    if(Health[xx, yy] == 0)
    {
      empty_cells = empty_cells + 1
      cells_coords[, empty_cells] = c(xx, yy)
    }
  }
  if(empty_cells > 0)
  {
    if (empty_cells >= 2)
      rd = sample(1:empty_cells, 2, replace = FALSE) # On ne veut pas 2 fois la même position du coup
    else
      rd = 1
    for(index in rd){
      new_cell_x = cells_coords[1, index]
      new_cell_y = cells_coords[2, index]
      Health[new_cell_x, new_cell_y] <<- original_health
      Type[new_cell_x, new_cell_y] <<- typeb
    }
  }
  
  # remove bacterie
  Kill_bacteria(x,y)
}


Comsumption<-function(x, y)
{
  # I G F L T C
  # 1 2 3 4 5 6
  
  type <- Type[x, y]
  # setup des goûts de la bactéries
  if     (type == 1) food_taste = 2:4
  else if(type == 2) food_taste = 1:4
  else               food_taste = 5:6
  
  xx = (x-1):(x+1)
  yy = (y-1):(y+1)
  # Join edge of the grid
  for(i in 1:3){
    xx[i] = ((xx[i] - 1)%%x_max) + 1
    yy[i] = ((yy[i] - 1)%%y_max) + 1
  }
  
  eat_value <- 0
  
  # loop over neighboring cells
  for(j in xx)
  {
    for(k in yy)
    {
      for(f in food_taste) # attention ordre de preference
      {
        if(Food[j, k, f] > 0)
        {
          Food[j, k, f] <<- Food[j, k, f] - 1
          eat_value = eat_value + 1
          
          if(type == 1)  # production Lactate par B
            Food[x, y, 5] <<- Food[x, y, 5] + 1
          
          break # limite à 1 par cellule
        }
      }
    }
  }
  Health[x, y] <<- Health[x, y] + eat_value
}


Adjust_lysozyme<-function(x, y)
{
  if(Type[x, y] == 2){ # Bactérie C
    Health[x, y] <<- Health[x, y] - lyzo_dmg
  }
}

Kill_bacteria<-function(x,y){
  Health[x, y] <<- 0
  Type[x, y] <<- 0
}

##################
## MAIN

Init_bacteria(n_bacts, r_bact)
Placing_food()

GIB = rep(NA, times=(timestamp+1))
num_bacteria = matrix(data = 0, ncol = (timestamp+1), nrow = 3)

B = sum(Type[,] == 1)
C = sum(Type[,] == 2)
D = sum(Type[,] == 3)

num_bacteria[,1] = c(B,C,D) 


for(t in 1:timestamp)
{
  if((t-1)%%plotfreq == 0){
    plot(Type, key=NULL, xlab=paste("time step = ",t-1), ylab='', col=color, axis.col=NULL, axis.row=NULL)
    legend(legend = c("B", "C", "D"), col = 1:3, x = "bottomright", fill=color[2:4])
  }
  
  # for(i in sample(1:x_max))
  for(i in (1:x_max))
  {
    # for(j in sample(1:y_max))
    for(j in (1:y_max))
    {
      if(Health[i, j] != 0)
        Comsumption(i, j)
      
      if(Health[i, j] >= DThreshold)
        Cell_Division(i, j)
      
      if(t > lyso_time)
        Adjust_lysozyme(i, j)
      
      if(Type[i, j] != 0){
        if(Health[i, j] > IThreshold)
          Health[i, j] = Health[i, j] - cost_of_living
        else
          Kill_bacteria(i,j)
        }
      if(Health[i, j] < 0)
        Health[i, j] <- 0
    }
  }
  
  if(t %% food_delay == 0)
    Add_food()
  
  B = sum(Type[,] == 1)
  C = sum(Type[,] == 2)
  D = sum(Type[,] == 3)
  num_bacteria[,t+1] = c(B,C,D) 
  
  # setup steady state bacteria number
  if(t == lyso_time){
    ss <- rowMeans(num_bacteria[,(t+2-steady_state_size):(t+1)])
    Bs <- ss[1]; Cs <- ss[2]; Ds <- ss[3]
  }
} # end main loop

ts = 1:(timestamp+1)

# construct GIB
for(t in ts){
  B = num_bacteria[1,t]
  C = num_bacteria[2,t]
  D = num_bacteria[3,t]
  GIB[t] = (B/Bs) / ((C+D)/(Cs+Ds)) - 1
}

########
### plots 
plot.default(ts, c(num_bacteria[1,]), type = "l", col = color[2], lwd = 1.4, ylim = c(40, max(num_bacteria)), xlab = "time step", ylab = "Nombre de bactéries")
lines(ts, num_bacteria[2,], type = "l", col = color[3], lwd = 1.4)
lines(ts, num_bacteria[3,], type = "l", col = color[4], lwd = 1.4)
if (lyzo_dmg != 0) 
  abline(v=lyso_time, col="darkorange", lty=2, lwd = 1.4)
grid()
legend(legend = c("B", "C", "D"), col = color[2:4], lty=1, lwd = 2, x = "topleft", bg="transparent")

plot.default(ts, GIB, type = "l", xlab = "time step", ylab = "GBI")
GIB_tmp = GIB; GIB_tmp[GIB<0] = 0
polygon(c(0,ts, timestamp+1), c(0,GIB_tmp, 0),col='#22CF22')
GIB_tmp = GIB; GIB_tmp[GIB>0] = 0
polygon(c(0,ts,timestamp+1), c(0,GIB_tmp,0),col='#CF2222')
grid()
if (lyzo_dmg != 0) 
  abline(v=lyso_time, col="darkorange", lty=2, lwd = 1.4)

