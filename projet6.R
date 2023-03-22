# ----- IMSVS ----- #

library('plot.matrix')

#This is is main body of the program 

## PARAM 
y_max = 50
x_max = 50
DThreshold = 50
IThreshold = 5
init_health = 20
n_bacts = 3

lyso_time = 50
division_factor = 0.5
init_food_value = 5
food_per_itt = 1 
food_delay = 3 # refill delay
lyzo_dmg = 5

timestamp = 500
steady_state_limit = timestamp / 10 # limite phase convergence, temps à partir duquel on recup le nbr moyen de bacterie pour contruire Bs Cs Ds 

color=c("#333333", "#BF75F9", "#35AEF7", "#CFCB3A")
plotfreq = 100

set.seed(1)
controle_test = TRUE

#!!!!!!!!!!!!!!!!!!!
    # TODO : note si steady-state = moyenne à la fin, comment on gere l'ajout de lysosyme vu que impacte le nb de bacterie
    # ---> on lance 2 experiences ? sans et avec lysosime, on recupere les Bs Cs Ds de la premiere
#!!!!!!!!!!!!!!!!!!!

## END PARAM

# Placing  bacteria
Init_bacteria<-function(n_bacts)
{
Health <<- matrix(data = 0, ncol =  x_max, nrow = y_max)
Type <<- matrix(data = 0, ncol =  x_max, nrow = y_max)

id = sample(1:(x_max*y_max), n_bacts, replace = FALSE)

  for (i in 1:n_bacts) {
    x = ((id[i]-1) %% y_max) + 1
    y = floor((id[i]-1) / y_max) + 1
    i = i + 1
    Health[x,y] <<- init_health
    
    # -> random
    # rd = runif(1)
    # if(rd < 0.33)      Type[x,y] <<- 1 # 'B'
    # else if(rd < 0.66) Type[x,y] <<- 2 # 'C'
    # else               Type[x,y] <<- 3 # 'D'
    
    # -> fixe
    if(i %% 3 == 0)      Type[x,y] <<- 1 # 'B'
    else if(i %% 3 == 1) Type[x,y] <<- 2 # 'C'
    else                 Type[x,y] <<- 3 # 'D'
  }
}

Init_bacteria(n_bacts)

# Placing the food
Placing_food<-function()
{
  Food <<- array(rep(init_food_value, x_max * y_max * 6), dim=c(x_max, y_max, 6))
}

Add_food<-function()
{
  Food <<- Food + array(rep(food_per_itt, x_max * y_max * 6), dim=c(x_max, y_max, 6))
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
      n_bacts <<- n_bacts + 1
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
  
  tmp <- 0
  total = 0
  for(j in (x-1):(x+1))
  {
    for(k in (y-1):(y+1))
    {
      # Join edge of the grid
      if(j < 1 || j > x_max)
        j = ((j - 1)%%x_max) + 1
      if(k < 1 || k > y_max)
        k = ((k - 1)%%y_max) + 1
      
      for(f in 1:6) # type food.
      {
        total = total + 1
        if(Food[j, k, f] > 0)
        {
          if(
            (type == 1 && (f >= 2 && f <= 4)) # 3 food
            || 
            (type == 2 && (f <= 4))           # 4 food
            ||
            (type == 3 && (f >= 5))           # 2 food
            )
          {
            Food[j, k, f] <<- Food[j, k, f] - 1
            tmp = tmp + 1
            # if(type == 1)
            #   Food[x, y, 5] <<- Food[x, y, 5] + 1 # production Lactate par B
          }
        }
      }
    }
  }
  Health[x, y] <<- Health[x, y] + tmp
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

Init_bacteria(n_bacts)
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
    plot(Type, key=NULL, xlab=paste("time = ",t-1), ylab='', col=color, axis.col=NULL, axis.row=NULL)
    legend(legend = c("B", "C", "D"), col = 1:3, x = "bottomright", fill=color[2:4])
  }

  for(i in sample(1:x_max))
  {
    for(j in sample(1:y_max))
    {
      if(Health[i, j] != 0)
        Comsumption(i, j)
      
      if(Health[i, j] > DThreshold)
        Cell_Division(i, j)
      
      Health[i, j] = Health[i, j] - IThreshold
        
      if(t > lyso_time && !controle_test)
        Adjust_lysozyme(i, j)
      
      if(Type[i, j] != 0 && Health[i, j] <= 0)
        Kill_bacteria(i,j)
      if(Health[i, j] <= 0)
        Health[i, j] <- 0
    }
  }
  
  if(t %% food_delay == 0)
    Add_food()
  
  B = sum(Type[,] == 1)
  C = sum(Type[,] == 2)
  D = sum(Type[,] == 3)
  num_bacteria[,t+1] = c(B,C,D) 
}

ts = 1:(timestamp+1)

# steady state : (moyenne en retirant la phase de convergence)
if (controle_test){
  ss <- rowSums(num_bacteria[,(timestamp+1-steady_state_limit):(timestamp+1)])/steady_state_limit
} else{
  GIB_control = GIB
}
Bs = ss[1]; Cs = ss[2]; Ds = ss[3]

# construct GIB
for(t in ts){
  B = num_bacteria[1,t]
  C = num_bacteria[2,t]
  D = num_bacteria[3,t]
  GIB[t] = (B/Bs) / ((C+D)/(Cs+Ds)) - 1
}
  
# plots 
plot.default(ts, c(num_bacteria[1,]), type = "l", col = color[2], lwd = 1.4, ylim = c(40, max(num_bacteria)), xlab = "time", ylab = "Nombre de bactéries")
lines(ts, num_bacteria[2,], type = "l", col = color[3], lwd = 1.4)
lines(ts, num_bacteria[3,], type = "l", col = color[4], lwd = 1.4)
if (!controle_test) 
  abline(v=lyso_time, col="darkorange", lty=2, lwd = 1.4)
abline(v=(timestamp+1-steady_state_limit), col="darkgreen", lty=2, lwd = 1.4) # limite de la phase de convergence pour le cacule de Bs,Cs,Ds
grid()
legend(legend = c("B", "C", "D"), col = color[2:4], lty=1, lwd = 2, x = "topleft", bg="transparent")

plot.default(ts, GIB, type = "l", xlab = "time")
GIB_tmp = GIB; GIB_tmp[GIB<0] = 0
polygon(c(0,ts, timestamp+1), c(0,GIB_tmp, 0),col='#22CF22')
GIB_tmp = GIB; GIB_tmp[GIB>0] = 0
polygon(c(0,ts,timestamp+1), c(0,GIB_tmp,0),col='#CF2222')
grid()
if (controle_test) 
  lines(ts, GIB_control, type = "l", xlab = "time", col="blue")
lines(ts, GIB, type = "l", xlab = "time")
lines(1:timestamp, rep(0,timestamp))
if (!controle_test) 
  abline(v=lyso_time, col="darkorange", lty=2, lwd = 1.4)

