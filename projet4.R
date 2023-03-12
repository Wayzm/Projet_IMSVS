# ----- IMSVS ----- #

#This is is main body of the program 

y_max = 50
x_max = 50
DThreshold = 50
IThreshold = 5
init_health = 20
t = 0
lyso_time = 10
division_factor = 0.25
init_food_value = 5
food_per_itt = 1
food_delay = 1
lyzo_dmg = 2
timestamp = 200
n_bacts = 30

# Placing  bacteria
Init_bacteria<-function(n_bacts)
{
coords <<- matrix(data = 0, ncol = n_bacts, nrow = 2)
Health <<- matrix(data = 0, ncol =  x_max, nrow = y_max)
Type <<- matrix(data = 0, ncol =  x_max, nrow = y_max)

id = sample(1:(x_max*y_max), n_bacts, replace = FALSE)

  for (i in 1:n_bacts) {
    x = ((id[i]-1) %% y_max) + 1
    y = floor((id[i]-1) / y_max) + 1
    coords[, i] <<- c(x, y)
    i = i + 1
    Health[x,y] <<- init_health
    rd = runif(1)
    if(rd < 0.33)      Type[x,y] <<- 1 # 'B'
    else if(rd < 0.66) Type[x,y] <<- 2 # 'C'
    else               Type[x,y] <<- 3 # 'D'
  }
}

Init_bacteria(n_bacts)
color = diag(Type[coords[1,],coords[2,]])
plot(coords[1,],coords[2,], xlim = c(1, x_max), ylim = c(1, y_max), col=color)

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
fin = 0
counter = 1
xs <- c(x - 1, x    , x + 1, x    , x - 1, x - 1, x + 1, x + 1)
ys <- c(y    , y + 1, y    , y - 1, y - 1, y + 1, y + 1, y - 1)

# We get the number of empty cells and their coords
empty_cells = 0 # Nombre de cellule vide 
cells_coords <- matrix(data = 0, ncol =  8, nrow = 2)

  while(counter <= 8)
  {
    xx = xs[counter]
    yy = ys[counter]
    if(xx > 0 && xx <= x_max && yy > 0 && yy <= y_max)
    {
        if(Health[xx, yy] == 0)
      {
          empty_cells = empty_cells + 1
          cells_coords[1, empty_cells] = xx
          cells_coords[2, empty_cells] = yy
        }
    }
    counter = counter + 1
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
      coords <<- cbind(coords, c(new_cell_x, new_cell_y))
    }
  }
    
  # remove bacterie // Bah du coup non?
  Kill_bacteria(x,y)
}

Comsumption<-function(x, y)
{
  # I G F L T C
  # 1 2 3 4 5 6
  
  tmp <- 0
  tmp = 0
  if(x > 1 && x < x_max && y > 1 && y < y_max)
    for(j in (x-1):(x+1))
    {
      for(k in (y-1):(y+1))
      {
      for(f in 1:6) # type food.
      {
        if(Food[j, k, f] > 0)
        {
          if(
            Type[x, y] == 1 && (f >= 2 && f <= 4)
            || 
            Type[x, y] == 2 && (f <= 4)
            ||
            Type[x, y] == 3 && (f >= 5)
            )
            Food[j, k, f] <<- Food[j, k, f] - 1
            tmp = tmp + 1
        }
      }
    }
  }
  Health[x, y] <<- Health[x, y] + tmp
  # print(tmp)
}


Adjust_lysozyme<-function(x, y)
{
  if(Type[x, y] == 2){ # Bactérie C
    Health[x, y] <<- Health[x, y] - lyzo_dmg
  }
}

Kill_bacteria<-function(x,y){
  index = 0
  
  a = which(coords[1,] == x)
  b = which(coords[2,] == y)
  index = intersect(a,b)
  
  Health[x, y] <<- 0
  Type[x, y] <<- 0
  
  if(index != 0){
    n_bacts <<- n_bacts - 1
    coords <<- coords[,-index]
  }
}

Init_bacteria(n_bacts)
Placing_food()


GIB = rep(NA, times=timestamp)
num_bacteria = matrix(data = 0, ncol = timestamp, nrow = 3)

for(t in 1:timestamp)
{
  color = diag(Type[coords[1,],coords[2,]])
  plot(coords[1,],coords[2,], xlim = c(1, x_max), ylim = c(1, y_max), col=color)
  title(paste("timestamp = ",t))
  legend(legend = c("B", "C", "D"), col = 1:3, x = "bottomright", lty=1)

  for(i in 1:x_max)
  {
    for(j in 1:y_max)
    {
      if(Health[i, j] != 0)
        Comsumption(i, j)
      
      if(Health[i, j] > DThreshold)
        Cell_Division(i, j)
      
      Health[i, j] = Health[i, j] - IThreshold
        
      if(t > lyso_time)
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
  Bs = sum(Health[Type[,] == 1])
  Cs = sum(Health[Type[,] == 2])
  Ds = sum(Health[Type[,] == 3])
  
  GIB[t] = (B/Bs) / ((C+D)/(Cs+Ds)) - 1
  num_bacteria[,t] = c(B,C,D) 
}

ts = 1:timestamp

plot(ts, c(num_bacteria[1,]), type = "o", col = 1, ylim = c(40, max(num_bacteria)))
lines(ts, num_bacteria[2,], type = "o", col = 2)
lines(ts, num_bacteria[3,], type = "o", col = 3)
abline(v=lyso_time)
legend(legend = c("B", "C", "D"), col = 1:3, lty=1, x = "topright")

up = GIB>=0
down = GIB<=0
half = (max(ts[up]) + min(ts[down]))/2
plot(ts, GIB, type = "o")
# polygon(x = c(min(ts[up]), ts[up],half),
#         y = c(0, GIB[up],0),
#         col = "red" )
# polygon(x = c(half,ts[down],max(ts[down])),
#         y = c(0,GIB[down],0),
#         col = "blue" )
lines(1:timestamp, rep(0,timestamp))
abline(v=lyso_time)

