# ----- IMSVS ----- #


#This is is main body of the program 

y_max = 50
x_max = 50
DThreshold = 50
IThreshold = 5
t = 0
lyso_time = 5

# Make it possible to be divided by 3
n_bacts = 30

# Placing  bacteria
Init_bacteria<-function(n_bacts)
{
  coords <<- matrix(data = 0, ncol = n_bacts, nrow = 2)
  Health <<- matrix(data = 0, ncol =  x_max, nrow = y_max)
  Type <<- matrix(data = 0, ncol =  x_max, nrow = y_max)
  
  ii = sample(1:x_max, n_bacts, replace = FALSE)
  jj = sample(1:y_max, n_bacts, replace = FALSE)
  for (i in 1:n_bacts) {
    x = ii[i]
    y = jj[i]
    coords[,i] <<- c(x,y)
    Health[x,y] <<- 20
    rd = runif(1)
    if(rd < 0.33)      Type[x,y] <<- 1 # 'B'
    else if(rd < 0.66) Type[x,y] <<- 2 # 'C'
    else               Type[x,y] <<- 3 # 'D'
  }
}


# Placing the food
Placing_food<-function(n_food)
{
  Food <<- matrix(data = 0, ncol =  x_max, nrow = y_max)
  Food_Value <<- matrix(data = 0, ncol =  x_max, nrow = y_max)
  
  for(i in 1: y_max)
  {
    for(j in 1:x_max)
    {
      vec = runif(1)
      Food_Value[i, j] <<- 20
      
      if(vec[1] < 0.16)
        Food[i, j] <<- "I"
      else if(vec[1] < 0.32)
        Food[i, j] <<- "G"
      else if(vec[1] < 0.48)
        Food[i, j] <<- "F"
      else if(vec[1] < 0.64)
        Food[i, j] <<- "L"
      else if(vec[1] < 0.80)
        Food[i, j] <<- "T"
      else
        Food[i, j] <<- "C"
    }
  }
}

Add_food<-function(n_food)
{
  xx = sample(1:x_max, n_food,replace = TRUE)
  yy = sample(1:y_max, n_food,replace = TRUE)
  
  for(ii in 1:n_food)
  {
    vec = runif(1)
    i = xx[ii]
    j = yy[ii]
    Food_Value[i, j] <<- 20
    
    if(vec[1] < 0.16)
      Food[i, j] <<- "I"
    else if(vec[1] < 0.32)
      Food[i, j] <<- "G"
    else if(vec[1] < 0.48)
      Food[i, j] <<- "F"
    else if(vec[1] < 0.64)
      Food[i, j] <<- "L"
    else if(vec[1] < 0.80)
      Food[i, j] <<- "T"
    else
      Food[i, j] <<- "C"
  }
}

# Cell Division
Cell_Division<-function(x, y, division_probability)
{
  p = division_probability
  original_health = Health[x, y] / 2
  new_cell_health = original_health * p
  typeb<-Type[x, y]
  fin = 0
  counter = 1
  xs <- c(x - 1, x    , x + 1, x    , x - 1, x - 1, x + 1, x + 1)
  ys <- c(y    , y + 1, y    , y - 1, y - 1, y + 1, y + 1, y - 1)
  
  while(fin < 2 && counter <= 8)
  {
    xx = xs[counter]
    yy = ys[counter]
    # In the paper's algorithm there was an OR statement
    # This made no sense as it allows overflow out of the grid if either c or d is good 
    if(xx > 0 && xx <= x_max && yy > 0 && yy <= y_max)
    {
      if(Health[xx, yy] == 0)
      {
        Health[xx, yy] <<- original_health
        Type[xx, yy] <<- typeb
        n_bacts <<- n_bacts + 1
        coords <<- cbind(coords, c(xx, yy))
        
        fin = fin + 1
      }
    }
    counter = counter + 1
  }
  # remove bacterie
  Kill_bacteria(x,y)
}

Comsumption<-function(x, y)
{
  r = 0
  for(j in x-1:x+1)
  {
    for(k in y-1:y+1)
    {
      # Ok here there's an issue since we could check the type before the for loops
      if((
        Type[x, y] == 1 && (Food[j, k] == "G" || Food[j, k] == "F" || Food[j, k] == "L")
        || 
        Type[x, y] == 2 && (Food[j, k] == "G" || Food[j, k] == "F" || Food[j, k] == "L" || Food[j, k] == "I")
        ||
        Type[x, y] == 3 && (Food[j, k] == "T" || Food[j, k] == "C")
        ) && Food_Value[j, k] != 0
      )
      {
        Food_Value[j, k] <<- Food_Value[j, k] - 1
        r = r + 1
        if(Food_Value[j, k] == 0)
          Food[j, k] <<- 0
      }
    }
  }
  Health[x, y] <<- Health[x, y] + r
}


Adjust_lysozyme<-function(x, y)
{
  if(Health[x, y] > IThreshold  && Type[x, y] == 2){ 
    Health[x, y] <<- Health[x, y] - 2
  }
}

Kill_bacteria<-function(x,y){
  index = 0
  
  a = which(coords[1,] == x)
  b = which(coords[2,] == y)
  index = intersect(a,b)
  
  if(index != 0){
    n_bacts <<- n_bacts - 1
    Health[x, y] <<- 0
    Type[x, y] <<- 0
    coords <<- coords[,-index]
  }
}


Init_bacteria(n_bacts)
Placing_food(x_max * y_max)


timestamp = 20
for(t in 1:timestamp)
{
  color = diag(Type[coords[1,],coords[2,]])
  plot(coords[1,],coords[2,], xlim = c(1, x_max), ylim = c(1, y_max), col=color)
  for(i in 1:x_max)
  {
    for(j in 1:y_max)
    {
      if(Health[i, j] != 0)
        Comsumption(i, j)
      
      if(Health[i, j] > DThreshold)
        Cell_Division(i, j, 1.0)
      
      if(Health[i, j] > IThreshold){
        Health[i, j] = Health[i, j] - IThreshold
        if(Health[i, j] <= 0)
          Kill_bacteria(i,j)
      }
        
      if(t > lyso_time)
        Adjust_lysozyme(i, j)
    }
  }
  Add_food(50)
}

# plot
color = diag(Type[coords[1,],coords[2,]])
plot(coords[1,],coords[2,], xlim = c(1, x_max), ylim = c(1, y_max), col=color)


