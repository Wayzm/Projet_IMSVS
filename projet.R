# ----- IMSVS ----- #

# Function to adjust when Lysozyme are added
# This function simulates the fact that adding lysozyme regulates C types bacteries
# n_id : id of a bacteria
# coords : a 2 * n_total matrix storing the X and Y coordinates of the cells
# Type : the grid with all cells' type stored
# IThreshold : threshold of health score at which cell death occurs
# Health : the grid with all cells' health value
Adjust_lysozyme<-function(x, y, Type, IThreshold, Health)
{
  typeb <- Type[x, y]
  if(Health[x][y] > IThreshold  & Type[x, y] == "C"){ Health[x, y] = Health[x, y] - 2} 
  return((H = Health))
}


# Consumption of nutrient from neighboring cells
# Food : the grid with all the food type
# Food_Value : the grid with all the food value
# n_id : id of a bacteria
# coords : a 2 * n_total matrix storing the X and Y coordinates of the cells
# Type : the grid with all cells' type stored
# Health : the grid with all cells' health value
Comsumption<-function(x, y, Food, Food_Value, Type, Health)
{
  r = 0
  for(j in x-1:x+1)
  {
    for(k in y-1:y+1)
    {
      # Ok here there's an issue since we could check the type before the for loops
      if(Type[x, y] == "B" || Type[x, y] == "C")
      {
        if(Food[j, k] == "G" || Food[j, k] == "F" || Food[j, k] == "L")
        {
          if(Food_Value[j, k] != 0)
          {
            Food_Value[j, k] = Food_Value[j, k] - 1
            r = r + 1
          }
        }
      }
      if(Type[x, y] == "C")
      {
        if(Food[j, k] == "I")
        {
          if(Food_Value[j, k] != 0)
          {
            Food_Value[j, k] = Food_Value[j, k] - 1
            r = r + 1
          }
        }
      }
      if(Type[x, y] == "D")
      {
        if(Food[j, k] == "T" || Food[j, k] == "C")
        {
          if(Food_Value[j, k] != 0)
          {
            Food_Value[j, k] = Food_Value[j, k] - 1 
            r = r + 1
          }
        }
      }
      Health[x, y] = r
    }
  }
  return(list(Value = Food_Value, H = Health))
}


# Division of the current cell if requirements are fulfilled
# n_id : id of a bacteria
# n_bacts : number of bacteria
# coords : a 2 * n_total matrix storing the X and Y coordinates of the cells
# Type : the grid with all cells' type stored
# Health : the grid with all cells' health value
Cell_Division<-function(x, y, Type, division_probability, Health, Xmax, Ymax)
{
  n_bcts = n_bcts + 1
  p = division_probability
  original_health = Health[x, y] / 2
  new_cell_health = original_health * p
  typeb<-Type[x, y]
  fin = 0
  counter = 0
  xs<-{x - 1; x; x + 1; x; x - 1; x + 1; x + 1}
  ys<-{y; y + 1; y; y - 1; y - 1; y + 1; y + 1; y - 1}
  c = 0
  d = 0
  index = 1
  
  while(fin < 2 & counter < 8)
  {
    c = xs[counter]
    d = ys[counter]
    # In the paper's algorithm there was an OR statement
    # This made no sense as it allows overflow out of the grid if either c or d is good 
    if(c >= 0 & c < Xmax & d >= 0 & d < Ymax)
    {
      if(Health[c, d] == 0)
      {
        Health[c, d] = original_health
        Type[c, d] = typeb
        
        #In the paper, it was only fin = 1 which would lead to making every cell the neighborhood a new cell
        fin = fin + 1
      }
    }
    if(index < 8)
    {
      index = index + 1
    } else {
      index = 0
    }
    counter = counter + 1
  }
  return(list(Ty = Type, H = Health))
}


#This is is main body of the program 

y_max = 50
x_max = 50
DThreshold = 50
IThreshold = 5
t = 0
lyso_time = 20

# Make it possible to be divided by 3
n_bacts = 30
coords = matrix(data = 0, ncol = n_bacts, nrow = 2)
Type = matrix(data = 0, ncol =  x_max, nrow = y_max)
Health = matrix(data = 0, ncol =  x_max, nrow = y_max)
Food = matrix(data = 0, ncol =  x_max, nrow = y_max)
Food_Value = matrix(data = 0, ncol =  x_max, nrow = y_max)

# Placing  bacteria
ii = sample(1:x_max, n_bacts, replace = FALSE)
jj = sample(1:y_max, n_bacts, replace = FALSE)
for(i in 1: n_bacts)
{
  x = ii[i]
  y = jj[i]
  coords[1, i] = x
  coords[2, i] = y
  Health[x, y] = 20
  if(i%%3 == 1)
  {
    Type[x, y] = "B"
  } else if(i%%3 == 2)
  {
    Type[x, y] = "C"
  } else {
    Type[x, y] = "D"
  }
}

# Placing the food
n_food = x_max * y_max

for(i in 1: y_max)
{
  for(j in 1:x_max)
  {
    vec = runif(1)
    Food_Value[x, y] = 20
    if(vec[1] < 0.16)
    {
      Food[i, j] = "I"
    } else if(vec[1] < 0.32 )
    {
      Food[i, j] = "G"
    } else if(vec[1] < 0.48){
      Food[i, j] = "F"
    } else if(vec[1] < 0.64)
    {
      Food[i, j] = "L"
    } else if(vec[1] < 0.80){
      Food[i, j] = "T"
    } else {
      Food[i, j] = "C"
    }
  }
}
timestamp = 100
while(t < 100)
{
  for(i in 1:x_max)
  {
    for(j in 1:y_max)
    {
      if(Health[i, j] != 0)
      {
        out = Comsumption(i, j, Food, Food_Value, Type, Health)
        Food_Value = out$Value
        Health = out$H
      }
      if(Health[i, j] > DThreshold)
      {
        out = Cell_Division(x, y, Type, 1.0, Health, x_max, y_max)
        Type = out$Ty
        Health = out$H
      }
      if(Health[i, j] > IThreshold)
      {
        Health[i, j] = Health[i, j] - IThreshold
      } else {
        Health[i, j] = 0
      }
      if(t > lyso_time){
        out = Adjust_lysozyme(i, j, Type, IThreshold, Health)
        Health = out$H
        
      }
    }
  }
  t = t + 1
}

